echo 32 4 16 10000 100 100 | rdmd -version=prep .\gen_doubly_good.d > in-16.txt
type in-16.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out16-0.txt
type in-16.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out16-1.txt
type in-16.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out16-2.txt
type in-16.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out16-3.txt
type in-16.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out16-4.txt
type in-16.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out16-5.txt
type in-16.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out16-6.txt
type in-16.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out16-7.txt
type in-16.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out16-8.txt
type in-16.txt | rdmd -version=exec .\gen_doubly_good.d | rdmd .\wcds.d > out16-9.txt

