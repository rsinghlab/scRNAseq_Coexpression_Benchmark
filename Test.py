import os

file_list = os.listdir("./Benchmark/")
file_list.remove("__init__.py")
for each in file_list:
    os.rename("./Benchmark//{}".format(each), "./Benchmark/(TODO){}".format(each))