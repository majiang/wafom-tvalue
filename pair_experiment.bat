echo %time%
echo 32 12 4 10000 100 100 | rdmd find_pair.d > pairs.csv
echo %time%
rdmd shift.d < pairs.csv > result.csv
echo %time%
merge.py
echo %time%
