npsimulation -D Example1.detector -E Example1.reaction -O Example1 -B run.mac
npanalysis --last-sim  -O Example1
root ShowResults.C
