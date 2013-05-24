cat test0.fa >tmp
for i in {1..19}
do
	cat test$i.fa >> tmp
done
