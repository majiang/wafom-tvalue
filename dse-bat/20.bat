echo 32 4 20 10000 100 100 | rdmd -version=prep .\gen_doubly_good.d > in-20.txt
type in-20.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out20-0.txt
type in-20.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out20-1.txt
type in-20.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out20-2.txt
type in-20.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out20-3.txt
type in-20.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out20-4.txt
type in-20.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out20-5.txt
type in-20.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out20-6.txt
type in-20.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out20-7.txt
type in-20.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out20-8.txt
type in-20.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out20-9.txt

