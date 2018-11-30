# gdal2tiles_tif
##基础环境
* 1.numpy
* 2.PIL
* 3.gdal>=2.1
* 4.pyhton

##使用说明：


##说明
* 1.-s :空间参考方式(EPSG:4326)
* 2.-p :投影方式 (mercator)
* 3.-z :层级 示例（'0-3'）
* 4.输入文件（tif）
* 5.输出文件夹

##示例
```
python gdal2tiles_tif.py -s EPSG:4326 -p mercator -z '0-3' /Volumes/pioneer/gdal_Demo/translateResult/Z_NAFP_C_BABJ_20180701091809_P_CLDAS_RT_ASI_0P0625_HOR-PRS-2018070109.tif /Volumes/pioneer/gdal_Demo/tilesResult/20181107_14_59

```