echo "in bash script"

python test.py
wait ${!}
echo "first wait is over"
./test2.sh
wait ${!}
echo "second wait is over"
python test.py
