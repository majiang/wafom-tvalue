echo 32 4 10 10000 100 100 | rdmd -version=prep .\gen_doubly_good.d > in-10.txt
type in-10.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out10-0.txt
type in-10.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out10-1.txt
type in-10.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out10-2.txt
type in-10.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out10-3.txt
type in-10.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out10-4.txt
type in-10.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out10-5.txt
type in-10.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out10-6.txt
type in-10.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out10-7.txt
type in-10.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out10-8.txt
type in-10.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out10-9.txt

