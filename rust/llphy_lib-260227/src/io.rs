//! Module defining various IO utilities
//! like [`read_file`] and [`read_file_into_global`].
use anyhow::Error;
use bumpalo::{Bump, collections::Vec};
use std::{
    fs::File,
    io::{ErrorKind, Read},
    path::Path,
};

/// Read the file at `path` into a byte vector
/// managed by the global allocator.
///
/// For non-persistent (and large!) files.
pub fn read_file_into_global(path: &Path) -> Result<std::vec::Vec<u8>, Error> {
    Ok(std::fs::read(path)?)
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
