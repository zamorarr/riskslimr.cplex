# Prepare your package for installation here.
# Use 'define()' to define configuration variables.
# Use 'configure_file()' to substitute configuration values.

parse_args <- function(args) {
  #m <- regexec("^--with-cplex-studio-dir=(.*)$", args[1])
  #r <- regmatches(args[1], m)[[1]]

  x <- strsplit(args, "=")
  varnames <- vapply(x, function(xi) xi[1], character(1))
  varnames <- gsub("^--", "", varnames)
  vars <- lapply(x, function(xi) trimws(xi[2]))
  setNames(vars, varnames)
}

# get configure args
args <- commandArgs(TRUE)
if (length(args) > 1 && args[1] == "configure") {
  args <- args[-1]
  args <- parse_args(args)
}

# get studio dir
studio_dir <- Sys.getenv("CPLEX_STUDIO_DIR")
has_env_var <- !identical(studio_dir, "")
has_configure_var <- "with-cplex-studio-dir" %in% names(args)

if (!has_env_var && !has_configure_var) {
  # R CMD INSTALL --preclean --no-multiarch --with-keep.source --configure-args="--with-cplex-studio-dir='/home/bobby/opt/ibm/ILOG/CPLEX_Studio_Community201'" riskslimr
  stop("must provide either a CPLEX_STUDIO_DIR environment variable or a --with-cplex-studio-dir=/path/to/dir argument to configure-args", call. = FALSE)
} else if (has_configure_var) {
  # no environment variable. look for configure dir in options
  studio_dir <- args[["with-cplex-studio-dir"]]
}

# determine studio directory
studio_dir <- normalizePath(studio_dir, mustWork = FALSE)
if (!dir.exists(studio_dir)) {
  stop(sprintf("directory %s does not exist\n", studio_dir), call. = FALSE)
}

cat(sprintf("using cplex studio dir: %s\n", studio_dir))

# generate directory params
cplex_dir <- file.path(studio_dir, "cplex")
concert_dir <- file.path(studio_dir, "concert")
arch <- list.files(file.path(cplex_dir, "bin"))

# generate include paths and linking paths
CPLEX_INCLUDE_PATH <- file.path(cplex_dir, "include")
CONCERT_INCLUDE_PATH <- file.path(concert_dir, "include")
CPLEX_LINK_PATH <- file.path(cplex_dir, "lib", arch, "static_pic")
CONCERT_LINK_PATH <- file.path(concert_dir, "lib", arch, "static_pic")

define(
  CXX_STD = "CXX11",
  #"PKG_CFLAGS" = "-m64 -fPIC -fno-strict-aliasing -fexceptions -DNDEBUG",
  "PKG_CXXFLAGS" = sprintf("$(SHLIB_OPENMP_CXXFLAGS) -I%s -I%s", CPLEX_INCLUDE_PATH, CONCERT_INCLUDE_PATH),
  #"PKG_LIBS" = sprintf("$(SHLIB_OPENMP_CXXFLAGS) $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS) -L%s -L %s -lconcert -lilocplex -lcplex -lm -lpthread", CPLEX_LINK_PATH, CONCERT_LINK_PATH)
  "PKG_LIBS" = sprintf("$(SHLIB_OPENMP_CXXFLAGS) $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS) -L%s -L %s -lconcert -lilocplex -lcplex -lm -lpthread", CPLEX_LINK_PATH, CONCERT_LINK_PATH)
)

#print(configure_database())


configure_file("src/Makevars.in")
#configure_file("src/Makevars.win.in")
