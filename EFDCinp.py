"""
code: utf-8
#################################
#     DEM to EFDC+ .inp data    #
#################################
#            Nie Qiyang         #
#            Hokkudai.          #
#            2023/1209          #
#################################
"""
import copy

def read_asc_file(file_path):
    data = []
    with open(file_path, 'r') as file:
        # 读取文件头信息
        ncols = int(file.readline().split()[1])
        nrows = int(file.readline().split()[1])
        xllcorner = float(file.readline().split()[1])
        yllcorner = float(file.readline().split()[1])
        cellsize = float(file.readline().split()[1])
        nodata_value = int(file.readline().split()[1])

        # 添加3行nodata_value到顶部
        for _ in range(3):
            data.append([nodata_value] * (ncols + 6))

        # 读取数据矩阵，自顶向下读，生成的二维列表是[y][x]方向，先存的x。
        for _ in range(nrows):
            row_data = list(map(int, file.readline().split()))
            # 在每一行的左右各添加3列-9999
            row_data = [nodata_value] * 3 + row_data + [nodata_value] * 3
            data.append(row_data)

        # 添加3行nodata_value到底部
        for _ in range(3):
            data.append([nodata_value] * (ncols + 6))

    # 更新行数和列数
    ncols = ncols+6
    nrows = nrows+6
    xllcorner = xllcorner -3.0*cellsize
    yllcorner = yllcorner -3.0*cellsize
    return ncols, nrows, xllcorner, yllcorner, cellsize, nodata_value, data

#生成cell.inp
def write_cell_inp(file_path, celllt_inp_path, ncols, nrows, data, nodata_value):
    with open(file_path, 'w') as file, open(celllt_inp_path, 'w') as file2:
        # 写入文件头信息
        file.write(f"C Cell.inp file, {ncols} columns and {nrows} rows\n")
        file.write("C Project:\n")
        file2.write(f"C Celllt.inp file, {ncols} columns and {nrows} rows\n")
        file2.write("C Project:\n")

        # 列取模
        count_10 = "".join(f"         {i}" for i in range(1, ncols//10+1))
        file.write("C    "+count_10+'\n')
        count_10 = "".join(f"         {i % 10}" for i in range(1, 11))
        file2.write("C    "+count_10+'\n')

        count_str = "".join(str(i % 10) for i in range(1, ncols + 1))
        file.write(f"C    {count_str}\n")
        file2.write(f"C    {count_str}\n")

        acell=0
        bcell=0
        data0 = copy.deepcopy(data)
        # 写入数据矩阵，处理无数据值为0，有数据值为5，边界值为9.
        for i in range(nrows):
            file.write(f"{nrows - i:3d}  ")
            file2.write(f"{nrows - i:3d}  ")
            for j in range(ncols):
                if data[i][j] == nodata_value :
                    if (i==0 or j==0 or i==nrows-1 or j==ncols-1):
                        file.write("0")
                        file2.write("0")
                        data0[i][j] = 0
                    elif ( i >0 and i<nrows-1 and j>0 and j<ncols-1 and \
                        data[i - 1][j - 1] == nodata_value and data[i - 1][j] == nodata_value and \
                        data[i - 1][j + 1] == nodata_value and data[i][j - 1] == nodata_value and \
                        data[i][j + 1] == nodata_value and data[i + 1][j - 1] == nodata_value and \
                        data[i + 1][j] == nodata_value and data[i + 1][j + 1] == nodata_value ):
                        file.write("0")
                        file2.write("0")
                        data0[i][j] = 0
                    else:
                        file.write("9")
                        file2.write("9")
                        data0[i][j] = 9
                        bcell=bcell+1
                else:
                    file.write("5")
                    file2.write("5")
                    data0[i][j] = 5
                    acell=acell+1
            file.write("\n")
            file2.write("\n")
    print(acell,bcell)
    return data0

#生成dxdy.inp, lxly.inp, depdat.inp
def write_dxdy_lxly_inp(dxdy_inp_path, lxly_inp_path, depdat_inp_path, gefdc__inp_path, ncols, nrows, data0, cellsize, depth, xllcorner, yllcorner):
    count = 0
    for row in data0:
        for value in row:
            if value == 5:
                count += 1
    with open(dxdy_inp_path, 'w') as file, open(lxly_inp_path, 'w') as file2, open(depdat_inp_path, 'w') as file3, open(gefdc__inp_path, 'w') as file4:
        # 写入文件头信息
        file.write(f"C DXDY.INP FILE, IN FREE FORMAT ACROSS COLUMNS for  {count} Active Cells\n")
        file.write("C Project:  EFDC+\n")
        file.write("C                                           BOTTOM                 Veg\n")
        file.write("C   I    J        DX        DY      DEPTH     ELEV     ZROUGH      TYPE\n")

        file2.write(f"C LXLY.INP FILE, IN FREE FORMAT ACROSS COLUMNS for  {count} Active Cells\n")
        file2.write("C Project:  EFDC+\n")
        file2.write("C                                                                             WIND\n")
        file2.write("C   I    J           X           Y       CUE       CVE       CUN      CVN  SHELTER\n")

        #dxdy内容
        CUE = 1.0
        CVE = 0.0
        CUN = 0.0
        CVN = 1.0
        for i in range(ncols):
            for j in range(nrows, 0, -1):
                if data0[j-1][i] == 5:
                    depthij = depth[j-1][i]
                    # 84.0m是洞爷湖测水深时基准平面高度
                    elev = 84.0 - depth[j-1][i]
                    file.write(f"{i + 1 :5d}" + f"{nrows - j + 1:5d}  " + f"{cellsize :10.3f}" + f"{cellsize :10.3f}" + f"{depthij :10.3f}" + f"{elev :10.3f}" + "  1.0000E-02\n")
                    file2.write(f"{i + 1 :5d}" + f"{nrows - j + 1:5d}  " + f"{xllcorner + (i + 0.5) * cellsize :12.1f}" + f"{yllcorner + (ncols - j + 0.5) * cellsize :12.1f}"\
                                 + f"{CUE :10.5f}" + f"{CVE :10.5f}" + f"{CUN :10.5f}" + f"{CVN :10.5f}" + "   1.00\n")
                    file3.write(f"{xllcorner + (i + 0.5) * cellsize :12.1f}" + f"{yllcorner + (ncols - j + 0.5) * cellsize :12.1f}"+ f"{depthij :10.3f}" + "\n")
                '''
                elif data0[j-1][i] == 9:
                    file4.write(f"{i + 1 :5d}" + f"{nrows - j + 1:5d}  " + f"{xllcorner + (i + 0.5) * cellsize :12.1f}" + f"{yllcorner + (ncols - j + 0.5) * cellsize :12.1f}"\
                                 "\n")
                '''
        for i in range(ncols):
            for j in range(nrows):
                if data0[j][i] == 9:
                    file4.write(f"{i + 1 :5d}" + f"{nrows - j :5d}  " + f"{xllcorner + (i + 0.5) * cellsize :12.1f}" + f"{yllcorner + (ncols - j - 0.5) * cellsize :12.1f}"\
                                 "\n")        
        #data0等二维列表都是先y后x。

#生成corners.inp
def write_corners_inp(corners_inp_path, tempb_inp_path,  ncols, nrows, data0, cellsize,  xllcorner, yllcorner):
    with open(corners_inp_path, 'w') as file, open(tempb_inp_path, 'w') as file2:
        # 写入文件头信息
        file.write(f"* CORNERS.INP - CELL CORNER COORDINATES USED FOR PLOTTING BY EFDC_EXPLORER (DSI)\n")
        file.write(f"*    I     J"+"          X              Y"*4+"\n")
        file2.write(f"* TEMPB.INP - CELL CORNER COORDINATES USED FOR PLOTTING BY EFDC_EXPLORER (DSI)\n")
        file2.write(f"*      I    J    TEMB TBEDTHK"+"\n")
        for j in range(nrows, 0, -1):
            for i in range(ncols):
                if data0[j-1][i] == 5:
                    file.write(f"{i + 1 :6d}" + f"{nrows - j + 1:6d}" \
                                + f"{xllcorner + i * cellsize :13.3f}" + f"{yllcorner + (ncols - j) * cellsize :13.3f}" \
                                + f"{xllcorner + i * cellsize :13.3f}" + f"{yllcorner + (ncols - j + 1) * cellsize :13.3f}" \
                                + f"{xllcorner + ( i + 1 ) * cellsize :13.3f}" + f"{yllcorner + (ncols - j + 1) * cellsize :13.3f}" \
                                + f"{xllcorner + ( i + 1 ) * cellsize :13.3f}" + f"{yllcorner + (ncols - j) * cellsize :13.3f}"+ "\n")
                    file2.write("     "+f"{i + 1 :4d}" + f"{nrows - j + 1:4d}"+"   25.00   10.00\n")
def write_show_inp(file_path):

    #i,j定义屏幕上输出显示的位置
    i=0
    j=0
    with open(file_path, 'w') as file:
        file.write("C  SHOW.INP - OUTPUT INFO - Version: 10.1\n")
        file.write("C  Project: \n")
        file.write("C \n")
        file.write("C   NSTYPE    NSHOWR    ISHOWC    JSHOWC    ISHPRT\n")
        file.write("C  \n")
        file.write("C   ZSSMIN    ZSSMAX   SSALMAX\n")
        file.write(f"         3      1000{i:10d}{j:10d}        50\n")
        file.write("  -500.000  1000.000  1000.000\n")

def main():
    path_main = r'C:\Users\13RNR-i9\Desktop\inp/'

    asc_file_path = path_main + 'depthUTM.asc'
    cell_inp_path = path_main + 'cell.inp'
    celllt_inp_path = path_main + 'celllt.inp'
    dxdy_inp_path = path_main + 'dxdy.inp'
    lxly_inp_path = path_main + 'lxly.inp'
    corners_inp_path = path_main + 'corners.inp'
    show_inp_path = path_main + 'show.inp'
    tempb_inp_path = path_main + 'tempb.inp'
    depdat_inp_path = path_main + 'depdat.inp'
    gefdc__inp_path = path_main + 'gefdc.inp'

    ncols, nrows, xllcorner, yllcorner, cellsize, nodata_value, data = read_asc_file(asc_file_path)
    depth = copy.deepcopy(data)
    data0 = write_cell_inp(cell_inp_path, celllt_inp_path, ncols, nrows, data, nodata_value)
    write_dxdy_lxly_inp(dxdy_inp_path, lxly_inp_path, depdat_inp_path, gefdc__inp_path,ncols, nrows, data0, cellsize, depth, xllcorner, yllcorner)
    write_corners_inp(corners_inp_path, tempb_inp_path,  ncols, nrows, data0, cellsize,  xllcorner, yllcorner)
    write_show_inp(show_inp_path)
    print(f"Conversion complete. Output file: {cell_inp_path}")


if __name__ == "__main__":
    main()
