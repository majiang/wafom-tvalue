echo 32 4 8 10000 100 100 | rdmd -version=prep .\gen_doubly_good.d > in-8.txt
type in-8.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out8-0.txt
type in-8.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out8-1.txt
type in-8.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out8-2.txt
type in-8.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out8-3.txt
type in-8.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out8-4.txt
type in-8.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out8-5.txt
type in-8.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out8-6.txt
type in-8.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out8-7.txt
type in-8.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out8-8.txt
type in-8.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out8-9.txt

