# gdal2tiles_tif
##运行环境
* python 版本：2.7 
## 基础环境
* numpy 版本：1.15.3
* PIL 版本：5.0.0
* gdal 版本：2.3.1

## 使用说明：


## 外部调用参数说明(带*为必须传的参数)
* -s :空间参考方式(EPSG:4326)
* -p :投影方式 (默认mercator)
* ***-z** :层级 示例（'0-3'）
* ***输入文件**（*.tif）
* ***输出文件夹**
* --export_nodata 输出重新设置nodata值(xxx)

## 示例
```
python gdal2tiles_tif.py -s EPSG:4326 -p mercator -z '0-3' /Volumes/pioneer/gdal_Demo/translateResult/Z_NAFP_C_BABJ_20180701091809_P_CLDAS_RT_ASI_0P0625_HOR-PRS-2018070109.tif /Volumes/pioneer/gdal_Demo/tilesResult/20181107_14_59

```