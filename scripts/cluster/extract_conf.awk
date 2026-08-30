# Print the single [[conf]] block of a profiling TOML whose [conf.backend]
# n_threads equals `want`.
#
# One PBS job runs one thread count, but configs/profiling.toml holds the whole
# sweep. Rather than keeping one file per thread count in sync by hand, each job
# carves its own single-config TOML out of the reference file, which stays the
# only place where grid size, Reynolds and niters are defined.
#
#   awk -v want=8 -f extract_conf.awk configs/profiling.toml > profiling_t8.toml

function flush_block() {
  if (buf != "" && nthreads == want)
    printf "%s", buf
}

BEGIN { want = want + 0; buf = ""; nthreads = -1; found = 0 }

/^[ \t]*\[\[conf\]\]/ {
  if (buf != "" && nthreads == want) found = 1
  flush_block()
  buf = $0 "\n"
  nthreads = -1
  next
}

{
  # Lines before the first [[conf]] (comments, stray keys) belong to no block
  # and are dropped: emitting them would prepend junk to the generated file.
  if (buf == "") next
  buf = buf $0 "\n"
  if ($0 ~ /^[ \t]*n_threads[ \t]*=/) {
    v = $0
    sub(/.*=[ \t]*/, "", v)
    nthreads = v + 0
  }
}

END {
  if (buf != "" && nthreads == want) found = 1
  flush_block()
  # A silent empty file would make properties_test fail much later with an
  # opaque "campo 'conf' mancante"; say what actually went wrong.
  if (!found) {
    printf "extract_conf.awk: no [[conf]] block with n_threads = %d\n", want > "/dev/stderr"
    exit 1
  }
}
