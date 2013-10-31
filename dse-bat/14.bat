echo 32 4 14 10000 100 100 | rdmd -version=prep .\gen_doubly_good.d > in-14.txt
type in-14.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out14-0.txt
type in-14.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out14-1.txt
type in-14.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out14-2.txt
type in-14.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out14-3.txt
type in-14.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out14-4.txt
type in-14.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out14-5.txt
type in-14.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out14-6.txt
type in-14.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out14-7.txt
type in-14.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out14-8.txt
type in-14.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out14-9.txt

