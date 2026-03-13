//! Some functions for this executable
//! that I couldn't bother organizing.

use std::{io::Write};
use bumpalo::{Bump, boxed::Box};

/// Type erase a writer, putting it into the given memory arena.
pub fn alloc_dyn_writer<'a>(writer: impl Write + 'static, arena: &'a Bump) -> Box<'a, dyn Write> {
    unsafe { Box::from_raw(arena.alloc(writer) as &mut dyn Write as *mut dyn Write) }
}