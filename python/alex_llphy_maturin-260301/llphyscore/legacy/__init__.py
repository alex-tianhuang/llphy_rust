import sys
import traceback
AMINOACIDS = "ACDEFGHIKLMNPQRSTVWY"

def load_named_pdb_statistics(dbpath: str, max_xmer: int):
    feature_tagABs = {
        "S2.SUMPI": ["srpipi", "lrpipi"],
        "S3.WATER.V2": ["Water", "Carbon"],
        # "ssL" also in original but seems not intentional
        # "ssH" in database but not in here... yikes...
        "S4.SSPRED": ["ssH", "ssE"],
        "S5.DISO": ["disL", "disS"],
        "S6.CHARGE.V2": ["srELEC", "lrELEC"],
        "S7.ELECHB.V2": ["sr_hb", "lr_hb"],
        "S8.CationPi.V2": ["srCATPI", "lrCATPI"],
        "S9.LARKS.V2": ["larkSIM", "larkFAR"],
    }
    pdb_statistics_names: "list[tuple[str, tuple[str, str]]]" = []
    pdb_statistics = []
    with Warner(10) as warner:
        for grid_name, (tag_sr, tag_lr) in feature_tagABs.items():
            grid_scorer = load_grid_scorer(dbpath, grid_name, tag_sr, tag_lr, max_xmer, warner)
            pdb_statistics_names.append((grid_name, (tag_sr, tag_lr)))
            pdb_statistics.append(grid_scorer)
    return pdb_statistics_names, pdb_statistics

def load_grid_scorer(dbpath: str, grid_name: str, tag_sr: str, tag_lr: str, max_xmer: int, warner: "Warner"):
    from argparse import Namespace

    grid_scorer = Namespace()
    subdir = f"{dbpath}/{grid_name}"
    filepath = f"{subdir}/PCON2.FREQS.wBOOTDEV"
    try:
        grid_scorer.pair_freq_db_x, grid_scorer.pair_freq_db_y = load_pair_freq_db(
            filepath, tag_sr, tag_lr, max_xmer, warner
        )
    except KeyboardInterrupt:
        raise
    except Exception as e:
        raise RuntimeError(
            "failed to load frequency pair DB @ {}".format(filepath)
        ) from e
    filepath = f"{subdir}/STEP6_PICKLES/SC_GRIDS.pickle4"
    try:
        grid_scorer.z_grid_db = load_z_grid_db(filepath, max_xmer)
    except KeyboardInterrupt:
        raise
    except Exception as e:
        raise RuntimeError("failed to load Z grid DB @ {}".format(filepath)) from e
    xmer_dir = f"{subdir}/STEP4_AVGnSDEVS"
    grid_scorer.avg_sdev_db = load_avg_sdev_db(xmer_dir, max_xmer)
    return grid_scorer

class Warner:
    limit: int
    warn_count: int
    entered: bool

    def __init__(self, limit: int) -> None:
        self.limit = limit
        self.warn_count = 0
    
    @staticmethod
    def no_limit():
        return Warner(sys.maxsize)

    def __enter__(self):
        return self

    def warn(self, msg):
        if self.warn_count < self.limit:
            print(msg, file=sys.stderr)
        self.warn_count += 1

    def __exit__(self, exc_type, exc, tb):
        if self.warn_count > self.limit:
            suppressed_count = self.warn_count - self.limit
            print(
                f"...\n{suppressed_count} other messages suppressed",
                file=sys.stderr,
            )


def load_pair_freq_db(filepath: str, tag_sr: str, tag_lr: str, max_xmer: int, warner: Warner):
    import math

    intermediate_x: "dict[int, dict[str, dict[str, list[tuple[float, float]]]]]" = {}
    intermediate_y: "dict[int, dict[str, dict[str, list[tuple[float, float]]]]]" = {}
    with open(filepath) as file:
        for line in file:
            def assert_true(b):
                if not b:
                    raise RuntimeError("failed to parse line: {}".format(line))

            def parse(f, s):
                try:
                    return f(s)
                except ValueError as e:
                    raise RuntimeError("failed to parse line: {}".format(line)) from e

            parts = (line := line.strip()).split()
            assert_true(len(parts) % 3 == 1)
            pair_key, *triplets = parts
            parts = pair_key.split("_")
            assert_true(len(parts) == 3)
            aa_x, separation, aa_y = parts
            separation = parse(int, separation)
            # assert_true(separation >= 1)
            if separation < 1:
                warner.warn(
                    f"ignoring line with non-positive {separation = }: {line}"
                )
                continue
            for ptype, freq, sdev in zip(triplets[::3], triplets[1::3], triplets[2::3]):
                assert_true((x := ptype.startswith("X_")) or ptype.startswith("Y_"))
                freq = parse(float, freq)
                sdev = parse(float, sdev)
                if x:
                    target = intermediate_x
                    first_aa = aa_x
                    second_aa = aa_y
                    tag = ptype.removeprefix("X_")
                else:
                    target = intermediate_y
                    first_aa = aa_y
                    second_aa = aa_x
                    tag = ptype.removeprefix("Y_")
                # assert_true((sr := tag == tag_sr) or tag == tag_lr)
                if not ((sr := tag == tag_sr) or tag == tag_lr):
                    available_tags = [tag_sr, tag_lr]
                    warner.warn(f"ignoring line with unrecognized {tag = } ({available_tags = }): {line}")
                if (entry1 := target.get(separation)) is None:
                    target[separation] = entry1 = {}
                if (entry2 := entry1.get(first_aa)) is None:
                    entry1[first_aa] = entry2 = {}
                if (entry3 := entry2.get(second_aa)) is None:
                    entry2[second_aa] = entry3 = [
                        (math.nan, math.nan),
                        (math.nan, math.nan),
                    ]
                if sr:
                    entry3[0] = (freq, sdev)
                else:
                    entry3[1] = (freq, sdev)
    pair_freq_db_x: "list[dict[str, dict[str, list[tuple[float, float]]]]]" = []
    pair_freq_db_y: "list[dict[str, dict[str, list[tuple[float, float]]]]]" = []
    for target, intermediate in (
        (pair_freq_db_x, intermediate_x),
        (pair_freq_db_y, intermediate_y),
    ):
        for i in range(max_xmer):
            separation = i + 1
            if (entry1 := intermediate.get(separation)) is None:
                x = len(pair_freq_db_x) == 0
                print(f"{x = }", intermediate.keys())
                raise RuntimeError(
                    "expected dense mapping for all possible residue separations, failed to find entry for separation {}".format(
                        separation
                    )
                )
            assert_is_full_aamap(entry1)
            for entry2 in entry1.values():
                assert_is_full_aamap(entry2)
                for entry3 in entry2.values():
                    for freq, sdev in entry3:
                        assert_true(not math.isnan(freq) and not math.isnan(sdev))
        target.extend(map(lambda tup: tup[1], sorted(intermediate.items())))
    return pair_freq_db_x, pair_freq_db_y


def load_z_grid_db(filepath: str, max_xmer: int):
    import pickle

    with open(filepath, "rb") as file:
        unpickled = pickle.load(file)

    def assert_true(b, expected, but):
        if not b:
            raise RuntimeError("expected {}, but {}".format(expected, but()))

    def narrow(any, ty):
        if not isinstance(any, ty):
            s = str(any)
            if len(s) > 80:
                s = s[:77] + "..."
            raise RuntimeError(
                "expected {}, but found {}: {}".format(
                    ty.__name__, type(any).__name__, s
                )
            )
        return any

    z_grid_db: "dict[str, list[list[tuple[float, list[tuple[float, tuple[int, int, int]]]]]]]" = {}
    unpickled_map = narrow(unpickled, dict)
    assert_is_full_aamap(unpickled_map)
    for key, entry1 in unpickled_map.items():
        aa = narrow(key, str)
        assert_true(
            aa in AMINOACIDS, "aminacid", lambda: "found non-aminoacid: {}".format(aa)
        )
        intermediate1: "list[tuple[int, list[tuple[float, list[tuple[float, tuple[int, int, int]]]]]]]" = []
        entry1 = narrow(entry1, dict)
        for key, entry2 in entry1.items():
            xmer = narrow(key, int)
            assert_true(
                xmer >= 1,
                "mapping of positive integer `xmer`",
                lambda: "found a non-positive integer: {}".format(xmer),
            )
            assert_true(
                xmer <= max_xmer,
                "dense mapping for all possible residue separations",
                lambda: "unexpectedly high `xmer`: {}".format(xmer),
            )
            intermediate2: "list[tuple[float, list[tuple[float, tuple[int, int, int]]]]]" = []
            intermediate1.append((xmer, intermediate2))
            entry2 = narrow(entry2, dict)
            for key, entry3 in entry2.items():
                sr_gridpoint = narrow(key, float)
                assert_true(
                    float.is_integer(sr_gridpoint * 2),
                    "gridpoints of half integer value",
                    lambda: "unexpected gridpoint: {}".format(sr_gridpoint),
                )
                intermediate3: "list[tuple[float, tuple[int, int, int]]]" = []
                intermediate2.append((xmer, intermediate3))
                entry3 = narrow(entry3, dict)
                for key, entry4 in entry3.items():
                    lr_gridpoint = narrow(key, float)
                    assert_true(
                        float.is_integer(lr_gridpoint * 2),
                        "gridpoints of half integer value",
                        lambda: "unexpected gridpoint: {}".format(lr_gridpoint),
                    )
                    total_sr_lr = narrow(entry4, list)
                    assert_true(
                        len(total_sr_lr) == 3,
                        "triplets of `(total, short-range, long-range)`",
                        lambda: "found unexpected list: {}".format(total_sr_lr),
                    )
                    total, sr, lr = total_sr_lr
                    total = narrow(total, int)
                    sr = narrow(sr, int)
                    lr = narrow(lr, int)
                    intermediate3.append((lr_gridpoint, (total, sr, lr)))
                intermediate3.sort()
            intermediate2.sort()
        intermediate1.sort()
        z_grid_db[aa] = [itm2 for _, itm2 in intermediate1]
    return z_grid_db


def load_avg_sdev_db(xmer_dir: str, max_xmer: int):
    avg_sdev_db: "dict[str, list[tuple[tuple[float, float], tuple[float, float]]]]" = {
        aa: [] for aa in AMINOACIDS
    }

    def load_avg_sdev_xmer_file(filepath):
        intermediate: "dict[str, tuple[tuple[float, float], tuple[float, float]]]" = {}
        with open(filepath) as file:
            for line in file:

                def assert_true(b):
                    if not b:
                        raise RuntimeError("failed to parse line: {}".format(line))

                def parse(f, s):
                    try:
                        return f(s)
                    except ValueError as e:
                        raise RuntimeError(
                            "failed to parse line: {}".format(line)
                        ) from e

                parts = (line := line.strip()).split()
                assert_true(len(parts) >= 7)
                aa, _, sr_avg, sr_sdev, _, lr_avg, lr_sdev, *_ = parts
                assert_true(aa in AMINOACIDS)
                sr_avg = parse(float, sr_avg)
                sr_sdev = parse(float, sr_sdev)
                lr_avg = parse(float, lr_avg)
                lr_sdev = parse(float, lr_sdev)
                intermediate[aa] = ((sr_avg, sr_sdev), (lr_avg, lr_sdev))
        assert_is_full_aamap(intermediate)
        return intermediate

    for xmer in range(1, max_xmer + 1):
        filepath = f"{xmer_dir}/PCON2.xmer{xmer}"
        try:
            intermediate = load_avg_sdev_xmer_file(filepath)
        except KeyboardInterrupt:
            raise
        except Exception as e:
            raise RuntimeError(
                "failed to load xmer avg/stdev statistics @ {}".format(filepath)
            ) from e
        for aa, target in avg_sdev_db.items():
            target.append(intermediate[aa])
    return avg_sdev_db


def assert_is_full_aamap(m):
    for key in m:
        if key not in AMINOACIDS:
            raise RuntimeError(
                "expected mapping of aminoacids, found non-aminoacid key: {}".format(
                    key
                )
            )
    if len(m) != 20:
        missing = str(aa for aa in AMINOACIDS if aa not in m)
        raise RuntimeError(
            "expected dense mapping of aminoacids, but missing entry/entries for `{}`".format(
                missing
            )
        )
