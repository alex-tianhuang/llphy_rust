//! Some utility functions for the project
//! that I couldn't bother organizing.
use {crate::datatypes::{FEATURE_NAMES, PAIR_NAMES_AND_FEATURE_NAMES}, anyhow::Error, bumpalo::{Bump, collections::Vec}, std::{fs::File, io::{ErrorKind, Read as _}, path::Path}};

/// Shorthand for a list search.
/// 
/// Since we have this:
/// ```
/// const PAIR_NAMES_AND_FEATURE_NAMES: /* ... */ = [
///     (pair_name_1, [feat_name_1a, feat_name_1b]),
///     (pair_name_2, [feat_name_2a, feat_name_2b]),
///     /* ... */
/// ]
/// ```
/// Find the pair_name and pair which contain the given `feat_name`.
pub fn lookup_f2p(feature_name: &str) -> Option<(&'static str, [&'static str; 2])> {
    PAIR_NAMES_AND_FEATURE_NAMES.iter().find(|(_, names)| names.contains(&feature_name)).copied()
}
/// Shorthand for a list search over an `(&str, T)` list.
/// 
/// Lookup by the name in the first slot of the tuple.
pub fn lookup_by_pair<'a, T>(data: &'a [(&'a str, T)], pair_name: &str) -> Option<&'a T> {
    data.iter().find(|(name, _)| *name == pair_name).map(|t| &t.1)
}
/// Shorthand for finding the index of a `feat_name` in `FEATURE_NAMES`.
pub fn lookup_findex(feat_name: &str) -> Option<usize> {
    FEATURE_NAMES.iter().enumerate().find(|(_, name)| **name == feat_name).map(|(i, _)| i)
}
/// Read the file at `path` into the given memory arena as a vector of bytes.
///
/// Dev note
/// --------
/// Most of this function is ripped from [`std::io::default_read_to_end`].
pub fn read_file<'a>(path: &Path, arena: &'a Bump) -> Result<Vec<'a, u8>, Error> {
    let mut file = File::open(path)?;
    let size_hint = file.metadata().map(|m| m.len() as usize).ok();
    let mut buf = Vec::with_capacity_in(size_hint.unwrap_or(0), arena);
    /// @tianh
    /// Maybe specific to the platform this was developed on,
    /// but I can't find the platform invariant export of this
    /// so I'm just hard coding it!
    const DEFAULT_BUF_SIZE: usize = 0x2000;
    /// @tianh
    /// `e.is_interrupted()` is private but this is what I think it should be equivalent to.
    fn is_interrupted(e: &std::io::Error) -> bool {
        matches!(e.kind(), ErrorKind::Interrupted)
    }
    /// @tianh
    /// `bumpalo::collections::AllocErr` does not implement `StdErr`
    fn map_alloc_err(_: bumpalo::collections::CollectionAllocErr) -> Error {
        Error::msg("memory allocation failed")
    }
    let start_cap = buf.capacity();
    // @std::io
    // Optionally limit the maximum bytes read on each iteration.
    // This adds an arbitrary fiddle factor to allow for more data than we expect.
    let mut max_read_size = size_hint
        .and_then(|s| {
            s.checked_add(1024)?
                .checked_next_multiple_of(DEFAULT_BUF_SIZE)
        })
        .unwrap_or(DEFAULT_BUF_SIZE);

    const PROBE_SIZE: usize = 32;

    fn small_probe_read(r: &mut File, buf: &mut Vec<u8>) -> std::io::Result<usize> {
        let mut probe = [0u8; PROBE_SIZE];

        loop {
            match r.read(&mut probe) {
                Ok(n) => {
                    // @std::io
                    // there is no way to recover from allocation failure here
                    // because the data has already been read.
                    buf.extend_from_slice(&probe[..n]);
                    return Ok(n);
                }
                // @tianh
                // Err(ref e) if e.is_interrupted() => continue,
                // Replacing commented out line above with line below
                Err(ref e) if is_interrupted(e) => continue,
                Err(e) => return Err(e),
            }
        }
    }

    // @std::io
    // avoid inflating empty/small vecs before we have determined that there's anything to read
    if (size_hint.is_none() || size_hint == Some(0)) && buf.capacity() - buf.len() < PROBE_SIZE {
        let read = small_probe_read(&mut file, &mut buf)?;

        if read == 0 {
            return Ok(buf);
        }
    }

    let mut consecutive_short_reads = 0;

    loop {
        if buf.len() == buf.capacity() && buf.capacity() == start_cap {
            // @std::io
            // The buffer might be an exact fit. Let's read into a probe buffer
            // and see if it returns `Ok(0)`. If so, we've avoided an
            // unnecessary doubling of the capacity. But if not, append the
            // probe buffer to the primary buffer and let its capacity grow.
            let read = small_probe_read(&mut file, &mut buf)?;

            if read == 0 {
                return Ok(buf);
            }
        }

        if buf.len() == buf.capacity() {
            // @std::io
            // buf is full, need more space
            buf.try_reserve(PROBE_SIZE).map_err(map_alloc_err)?;
        }

        // @tianh
        // The proper thing to do is to use BorrowedBuf
        // However, I can't be bothered and I have been convinced by chatGPT
        // that most `File::read` implementations will not depend
        // on reading the contents of the buf that was passed in.
        // So here are some spooky uninitialized bytes!
        let mut spare =
            unsafe { std::slice::from_raw_parts_mut(buf.as_mut_ptr(), buf.capacity() - buf.len()) };
        let buf_len = std::cmp::min(spare.len(), max_read_size);
        spare = &mut spare[..buf_len];
        let bytes_read = file.read(spare)?;

        // @tianh
        // SAFETY: yep
        unsafe {
            let new_len = bytes_read + buf.len();
            buf.set_len(new_len);
        }

        if bytes_read == 0 {
            return Ok(buf);
        }

        if bytes_read < buf_len {
            consecutive_short_reads += 1;
        } else {
            consecutive_short_reads = 0;
        }

        // @std::io
        // Use heuristics to determine the max read size if no initial size hint was provided
        if size_hint.is_none() {
            // @std::io
            // The reader is returning short reads but it doesn't call ensure_init().
            // In that case we no longer need to restrict read sizes to avoid
            // initialization costs.
            // When reading from disk we usually don't get any short reads except at EOF.
            // So we wait for at least 2 short reads before uncapping the read buffer;
            // this helps with the Windows issue.
            if consecutive_short_reads > 1 {
                max_read_size = usize::MAX;
            }

            // @std::io
            // we have passed a larger buffer than previously and the
            // reader still hasn't returned a short read
            if buf_len >= max_read_size && bytes_read == buf_len {
                max_read_size = max_read_size.saturating_mul(2);
            }
        }
    }
}


/// Leak a [`Vec`] managed by an arena into a slice
/// that lives for as long as the arena does not reset.
pub fn leak_vec<'a, T>(buf: bumpalo::collections::Vec<'a, T>) -> &'a mut [T] {
    let mut buf = std::mem::ManuallyDrop::new(buf);
    // SAFETY: the arena must de-allocate before this vec does.
    //         Also, do to the use of `ManuallyDrop`, these bytes
    //         are not de-allocated and reused by the arena.
    unsafe { std::mem::transmute::<&mut [T], &'a mut [T]>(buf.as_mut_slice()) }
}

/// Boilerplate that looks like `#[serde(from)]`.
macro_rules! derive_borsh_de_from {
    ($my_ty:ty as $std_ty:ty, $from_impl:expr) => {
        impl borsh::BorshDeserialize for $my_ty {
            fn deserialize(buf: &mut &[u8]) -> std::io::Result<Self> {
                <$std_ty as borsh::BorshDeserialize>::deserialize(buf).map($from_impl)
            }
            fn deserialize_reader<R: std::io::Read>(reader: &mut R) -> std::io::Result<Self> {
                <$std_ty as borsh::BorshDeserialize>::deserialize_reader(reader).map($from_impl)
            }
        }
    };
}

/// Boilerplate that looks like `#[serde(into)]`.
macro_rules! derive_borsh_se_into {
    ($my_ty:ty as $std_ty:ty, |$this:ident| $to_impl:expr) => {
        impl borsh::BorshSerialize for $my_ty {
            fn serialize<W: std::io::Write>(&self, writer: &mut W) -> std::io::Result<()> {
                let $this = self;
                let intermediate: $std_ty = $to_impl;
                intermediate.serialize(writer)
            }
        }
    };
}
pub(crate) use {derive_borsh_de_from, derive_borsh_se_into};
