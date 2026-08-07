#' @keywords internal
.onLoad <- function(libname, pkgname) {
   # required so S7 methods registered on external generics (e.g. print)
   # are activated correctly when the package is loaded/lazy-loaded.
   S7::methods_register();
}
