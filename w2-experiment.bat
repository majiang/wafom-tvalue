dmd refer_sobol.d bf.d lib\pointset.d lib\sobol.d lib\graycode.d ui\input.d lib\wafom.d -J.
dmd beat_sobol.d bf.d lib\pointset.d lib\sobol.d lib\graycode.d ui\input.d lib\wafom.d
dmd rmse.d lib\wafom.d lib\pointset.d lib\graycode.d lib\sobol.d lib\integral.d lib\integration_error.d ui\input.d tf.d

refer_sobol 4 20 4 > sobol-w2.csv
beat_sobol 4 12 4 > low-w2-small.csv
beat_sobol 13 16 4 > low-w2-medium.csv
beat_sobol 17 18 4 > low-w2-large.csv
beat_sobol 19 19 4 > low-w2-larger.csv
beat_sobol 20 20 4 > low-w2-largest.csv
rmse < sobol-w2.csv > sobol-err.csv
rmse < low-w2-small.csv > low-err-small.csv
rmse < low-w2-medium.csv > low-err-medium.csv
rmse < low-w2-large.csv > low-err-large.csv
rmse < low-w2-larger.csv > low-err-larger.csv
rmse < low-w2-largest.csv > low-err-largest.csv
