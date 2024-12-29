import pandas as pd
import matplotlib.pyplot as plt
import os
import shutil

#make a folder called 'wfun'
def mkdir(path):
    path = path.strip()
    path = path.rstrip("\\")
    isExists = os.path.exists(path)
    if not isExists:
        os.makedirs(path)
        return True
    else:
        return False

def plot_cluster(file_path, sheet_name,day):
    # 读取Excel文件中特定的工作表
    if day=='Day1':
        df = pd.read_excel(file_path, sheet_name=sheet_name, usecols="B,C", skiprows=1, skipfooter=1)
    elif day=='Day2':
        df = pd.read_excel(file_path, sheet_name=sheet_name, usecols="E,F", skiprows=1, skipfooter=1)
    elif day=='Day3':
        df = pd.read_excel(file_path, sheet_name=sheet_name, usecols="H,I", skiprows=1, skipfooter=1)
    elif day=='Special':
        df = pd.read_excel(file_path, sheet_name=sheet_name, usecols="K,L", skiprows=1, skipfooter=1)

    # 删除空白行
    df.dropna(how='all', inplace=True)

    # 将第一行非空白数据作为标题
    df.columns = df.iloc[1]
    df = df[2:]

    # 重置索引
    df.reset_index(drop=True, inplace=True)

    name_option_dict = dict(zip(df['name'], df['option']))
    return name_option_dict

# Copy file to another folder based on the name
def cp_file(name, option):
    print(name)
    if option == 'K':
        file_name = f'norm_wfun_{name}.fits'
    elif option == 'M':
        file_name = f'norm_M_wfun_{name}.fits'
    else:
        print(f"Unknown option for {name}")
        return

    source_file = os.path.join(source_folder, file_name)
    if os.path.exists(source_file):
        shutil.copy(source_file, target_folder)
        print(f"Copied: {file_name}")
    else:
        print(f"File not found: {file_name}")

def delete(target_folder):
    # Define the directory path to be deleted

    # Ensure the directory exists before attempting to delete
    if os.path.exists(target_folder):
        shutil.rmtree(target_folder)

# 使用示例
day='Day1'
day='Special'
file_path = r'C:\Users\ZY\Documents\github\233boy\Dr.-Yep-2024-summer-research\Goodman_Observing_Log.xlsx'
sheet_name = 'selection'

if day=='Special':
    source_folder = r'C:\Users\ZY\Documents\github\233boy\Dr.-Yep-2024-summer-research\Special\special_efiles\wfun\norm'
    target_folder =  r'C:\Users\ZY\Documents\github\233boy\Dr.-Yep-2024-summer-research\Special\special_efiles\wfun\norm_final'
else:
    source_folder = f"C:\\Users\\ZY\\Documents\\github\\233boy\\Dr.-Yep-2024-summer-research\\{day}\\RED\\wfun\\norm"
    target_folder = f"C:\\Users\\ZY\\Documents\\github\\233boy\\Dr.-Yep-2024-summer-research\\{day}\\RED\\wfun\\norm_final"

delete(target_folder)
mkdir(target_folder)
dic=plot_cluster(file_path, sheet_name, day)
print(dic)
# name='HD32450'
# option='M'
# cp_file(name,option)
for name,option in plot_cluster(file_path, sheet_name, day).items():
    cp_file(name,option)
