echo 32 4 12 10000 100 100 | rdmd -version=prep .\gen_doubly_good.d > in-12.txt
type in-12.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out12-0.txt
type in-12.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out12-1.txt
type in-12.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out12-2.txt
type in-12.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out12-3.txt
type in-12.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out12-4.txt
type in-12.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out12-5.txt
type in-12.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out12-6.txt
type in-12.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out12-7.txt
type in-12.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out12-8.txt
type in-12.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out12-9.txt

