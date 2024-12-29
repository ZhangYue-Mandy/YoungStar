
import pandas as pd
import matplotlib.pyplot as plt
from PIL import Image, ImageDraw, ImageFont
import os
import shutil
from astropy.io import fits
import re

def openfits(fitsfile_fp):
    fitsfile = fits.open(fitsfile_fp)
    data = fitsfile[0].data
    #head=fitsfile[0].header
    w=data[0]
    f=data[1]
    fitsfile.close()
    return w,f


def plot_cluster(file_path, sheet_name):
    # 读取Excel文件中特定的工作表
    df = pd.read_excel(file_path, sheet_name=sheet_name, usecols="E,F", skiprows=0, skipfooter=0)
    # 删除空白行
    df.dropna(how='all', inplace=True)

    # 将第一行非空白数据作为标题
    df.columns = df.iloc[0]
    df = df[1:]
    df.reset_index(drop=True, inplace=True)
    # Debugging: Print the DataFrame columns
    print("DataFrame columns:", df.columns)

    # Ensure column names are correct
    df.columns = [col.strip() for col in df.columns]

    name_option = zip(df['object'], df['SpT'])
    return name_option




#the cut range
# 6275 - 6315
# 6470 - 6530
# 6570 - 6580
def plot_skip(w, f, name):
    # Cut ranges
    cut_ranges = [(6275, 6315), (6470, 6530), (6570, 6580)]

    # Filter the data
    filtered_w = []
    filtered_f = []

    for wi, fi in zip(w, f):
        if not any(start <= wi <= end for start, end in cut_ranges):
            filtered_w.append(wi)
            filtered_f.append(fi)

    # Plot the filtered data
    plt.plot(w,f,color='grey',alpha=0.5)
    plt.plot(filtered_w, filtered_f)
    plt.xlabel('w')
    plt.ylabel('f')
    plt.ylim(0.3,1.5)
    plt.annotate(name, xy=(0.01, 1.02), xycoords='axes fraction', ha='left', fontsize=12)
    plt.savefig(output_path+'\\'+starname+'.png')
    plt.show()
    

def extract_filename(file_path):
    return file_path[file_path.rfind('\\') + 1:]

def extract_target(filename):
    if filename.endswith('.fits'):
        return filename[0:-len('.fits')]
def sort(elements, filenames):
    # 清理元素和文件名中的空格
    clean_elements = [element.strip() for element in elements]
    clean_filenames = {filename: filename for filename in filenames}

    # 创建一个映射元素到文件名的字典
    element_to_filename = {element: filename for filename in clean_filenames for element in clean_elements if element in filename}

    # 按照元素列表的顺序对文件名进行排序
    sorted_filenames = [element_to_filename[element] for element in clean_elements]

    return sorted_filenames
def combine_images_vertical(image_paths, output_path):
    images = [Image.open(path) for path in image_paths]
    width, height = images[0].size
    
    # Calculate the size of the combined image
    combined_width = width
    combined_height = height * len(images)
    
    combined_image = Image.new('RGB', (combined_width, combined_height))
    
    for index, image in enumerate(images):
        y_offset = index * height
        combined_image.paste(image, (0, y_offset))
    
    combined_image.save(output_path+'\\combined.jpg')
    
def combine_images_grid(image_paths, output_path, grid_size):
    images = [Image.open(path) for path in image_paths]
    width, height = images[0].size
    
    # Calculate the size of the combined image
    combined_width = width * grid_size[0]
    combined_height = height * grid_size[1]
    
    combined_image = Image.new('RGB', (combined_width, combined_height))
    
    for index, image in enumerate(images):
        x_offset = (index % grid_size[0]) * width
        y_offset = (index // grid_size[0]) * height
        combined_image.paste(image, (x_offset, y_offset))
    
    combined_image.save(output_path)


def stack(sorted_image_paths, standard_img):
    for i in range(len(sorted_image_paths)):
        adapted_paths=[]
        for _ in range(2):
            adapted_paths.append(standard_img[_])
        for _ in range(2):
            adapted_paths.append(sorted_image_paths[i])
        for _ in range(2):
            adapted_paths.append(standard_img[_+2])
        for _ in range(2):
            adapted_paths.append(sorted_image_paths[i])
        adapted_paths.append(standard_img[_+3])
        adapted_paths.append(sorted_image_paths[i])
        adapted_paths.append(sorted_image_paths[i])
        print(sorted_image_paths[i])
        combine_images_grid(adapted_paths, output_path+'\\'+extract_filename(extract_target(sorted_image_paths[i]))+'.jpg',[2,6])
    return adapted_paths

def scan_file(source_folder, extension=".png"):
    file_list = os.listdir(source_folder)
    specific_files = [os.path.join(source_folder, file) for file in file_list if file.lower().endswith(extension.lower())]
    return specific_files

file_path=r'C:\Users\ZY\Documents\github\233boy\Dr.-Yep-2024-summer-research\Standard_star.xlsx'
output_path=r'C:\Users\ZY\Documents\github\233boy\Dr.-Yep-2024-summer-research\CHIRONstds'


# #read the excel to get the name and spectral type

# sheet_name='Sheet3'
# dic=list(plot_cluster(file_path, sheet_name))
# files=scan_file(output_path)
# elements=[dic[i][0].strip() for i in range(len(dic))]
# sorted_filenames = sort(elements, files)

files=scan_file(output_path)
#plot the spectrum

# for i in files:
#     print(i)
#     w,f=openfits(i)
#     name=extract_target(extract_filename(i))
#     starname=extract_target(extract_filename(i))
#     plot_skip(w, f, name)

for i in files:
    combine_images_vertical(files, output_path)

# files=scan_file(output_path)
# #plot the spectrum

# for i in files:
#     print(i)
#     w,f=openfits(i)
#     name=extract_target(i)
#     starname=extract_target(extract_filename(i))
#     plot_skip(w, f, name)

# image_paths= scan_file(r'C:\Users\ZY\Documents\github\233boy\Dr.-Yep-2024-summer-research\norm_final0',extension='.png')
# standard_img=scan_file(r'C:\Users\ZY\Documents\github\233boy\Dr.-Yep-2024-summer-research\standard',extension='.png')

# # sheet_name1='Sheet2'
# # dic1=list(plot_cluster(file_path, sheet_name1))
# # elements1=[dic1[i][0].strip() for i in range(len(dic1))]
# # sorted_standard_img = sort(elements1, standard_img)

# # sorted_image_paths= sort(elements, image_paths)

# stack(image_paths, standard_img) 


