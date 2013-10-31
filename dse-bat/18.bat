echo 32 4 18 10000 100 100 | rdmd -version=prep .\gen_doubly_good.d > in-18.txt
type in-18.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out18-0.txt
type in-18.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out18-1.txt
type in-18.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out18-2.txt
type in-18.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out18-3.txt
type in-18.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out18-4.txt
type in-18.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out18-5.txt
type in-18.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out18-6.txt
type in-18.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out18-7.txt
type in-18.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out18-8.txt
type in-18.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out18-9.txt

