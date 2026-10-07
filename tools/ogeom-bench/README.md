# ogeom-bench

Wall-clock benchmarks over the kernel's hot paths.

```sh
cargo run --release -p ogeom-bench                     # table on stderr, baseline JSON on stdout
cargo run --release -p ogeom-bench -- --check tools/ogeom-bench/baseline.json
cargo run --release -p ogeom-bench -- --threads 4 --filter boolean
```

Each benchmark runs once to warm up, then repeats until at least 50 ms of
samples are timed (at least 3 samples, at most about 3 s of wall clock).
The table gives the minimum, the median, the median absolute deviation
(MAD) and the sample count. The whole run takes about a minute at one
thread.

## Ratios

Every run also times a fixed single-threaded arithmetic spin. A benchmark's
ratio is its minimum time over the spin's minimum, so it carries between
machines where milliseconds do not. The minimum is compared because load on
a shared machine only adds time. Even so, two runs on an idle machine
differ by up to about 10 percent, so a drift smaller than that is noise.

## Thread counts

The kernel's parallel stages use `ogeom_core::parallel::threads()`: the
count given by `--threads N`, else the `OGEOM_THREADS` environment
variable, else the machine's parallelism.

The spin is single-threaded, so a ratio only means the same thing at the
same thread count. `baseline.json` records the count it was taken at
(`"threads": 1`) and `--check` runs at that count unless `--threads` says
otherwise, in which case it warns that the ratios do not compare. At one
thread the ratios compare across machines whatever their core count. To
measure a parallel stage, compare paired runs at the same `--threads N` on
one machine.

## Renewing the baseline

Run at the lowest load you can find, three times, and keep each
benchmark's smallest ratio:

```sh
for r in 1 2 3; do cargo run --release -p ogeom-bench -- --threads 1 > run$r.json; done
```

The current baseline was taken at one thread on a 20-thread machine with a
load average of about 2.
