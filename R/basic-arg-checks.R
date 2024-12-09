
lgl_check <- function(x, n) {
    if (!(is.logical(x) && length(x) == 1)) {
        stop(paste("\nERROR:", n, "is not a single logical.\n"))
    }
}

uint_check <- function(x, n, .max = NULL, .min = NULL) {
    if (!(is.numeric(x) && length(x) == 1 && x %% 1 == 0 && x >= 0)) {
        stop(paste("\nERROR:", n, "cannot be properly cast as an",
                   "unsigned integer.\n"))
    }
    if (!is.null(.min) && x < .min) {
        stop(paste0("\n", n, " is below the minimum allowed value (",
                    .min, ").\n"))
    }
    if (!is.null(.max) && x > .max) {
        stop(paste0("\n", n, " is above the maximum allowed value (",
                    .max, ").\n"))
    }
}
uint_vec_check <- function(x, n, .max = NULL, .min = NULL) {
    if (!(is.numeric(x) && all(x %% 1 == 0) && all(x >= 0))) {
        stop(paste("\nERROR:", n, "cannot be properly cast as an",
                   "unsigned integer vector.\n"))
    }
    if (!is.null(.min) && any(x < .min)) {
        stop(paste0("\n", n, " contains values below the minimum ",
                    "allowed (", .min, ").\n"))
    }
    if (!is.null(.max) && any(x > .max)) {
        stop(paste0("\n", n, " contains values above the maximum ",
                    "allowed (", .max, ").\n"))
    }
}
dbl_check <- function(x, n, .max = NULL, .min = NULL) {
    if (!(is.numeric(x) && length(x) == 1)) {
        stop(paste("\nERROR:", n, "cannot be properly cast as a",
                   "double.\n"))
    }
    if (!is.null(.min) && x < .min) {
        stop(paste0("\n", n, " is below the minimum allowed value (",
                    .min, ").\n"))
    }
    if (!is.null(.max) && x > .max) {
        stop(paste0("\n", n, " is above the maximum allowed value (",
                    .max, ").\n"))
    }
}
dbl_vec_check <- function(x, n, .max = NULL, .min = NULL) {
    if (!(is.numeric(x) && is.null(dim(x)))) {
        stop(paste("\nERROR:", n, "cannot be properly cast as a",
                   "numeric vector.\n"))
    }
    if (!is.null(.min) && any(x < .min)) {
        stop(paste0("\n", n, " contains values below the minimum ",
                    "allowed (", .min, ").\n"))
    }
    if (!is.null(.max) && any(x > .max)) {
        stop(paste0("\n", n, " contains values above the maximum ",
                    "allowed (", .max, ").\n"))
    }
}
dbl_mat_check <- function(x, n, .max = NULL, .min = NULL) {
    if (!(is.numeric(x) && inherits(x, "matrix"))) {
        stop(paste("\nERROR:", n, "cannot be properly cast as a",
                   "numeric matrix.\n"))
    }
    if (!is.null(.min) && any(x < .min)) {
        stop(paste0("\n", n, " contains values below the minimum ",
                    "allowed (", .min, ").\n"))
    }
    if (!is.null(.max) && any(x > .max)) {
        stop(paste0("\n", n, " contains values above the maximum ",
                    "allowed (", .max, ").\n"))
    }
}
cube_list_check <- function(x, n) {
    if (!(inherits(x, "list") &&
          all(sapply(x, inherits, what = "array")) &&
          all(sapply(x, function(y) length(dim(y)) == 3)))) {
        stop(paste("\nERROR:", n, "cannot be properly cast as a",
                   "list of cubes.\n"))
    }
}

