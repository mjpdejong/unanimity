
Test doing 0-pass ccs on a subread from the lexogen-SIRV dataset

  $ $__PBTEST_CCS_EXE --logLevel=DEBUG --minPasses=0 --minPredictedAccuracy=0.85 $TESTDIR/../data/0passes.bam test.fq
  >|> \d{8} \d{2}:\d{2}:\d{2}\.\d{3} -|- DEBUG      -|- main -|- [0-9,a-f,x]+|| -|- Found consensus models for: (P6-C4, S/P1-C1) (re)
  >|> \d{8} \d{2}:\d{2}:\d{2}\.\d{3} -|- DEBUG      -|- main -|- [0-9,a-f,x]+|| -|- Using consensus models for: (P6-C4) (re)
  $ cat test.fq
  cat: test.fq: No such file or directory
  [1]
