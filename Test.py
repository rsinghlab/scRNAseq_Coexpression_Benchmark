import os

file_list = os.listdir("./Plotting/")
file_list.remove("__init__.py")
for each in file_list:
    os.rename("./Plotting//{}".format(each), "./Plotting/(TODO){}".format(each))