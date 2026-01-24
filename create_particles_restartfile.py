#!/usr/bin/env python
# coding: utf-8
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import netCDF4 as nc

def createParticlesNC(fnam,ids,ines,jnes,lons,lats,slons,slats,us,vs,depths,kfixed,days):
    f=nc.Dataset(fnam,'w',format='NETCDF3_64BIT_OFFSET',clobber=True)
    idim=f.createDimension('i',None)
    ivv=f.createVariable('i',np.float64)
    iv=f.createVariable('drifter_num',np.int32,('i',))
    inev=f.createVariable('ine',np.int32,('i',))
    jnev=f.createVariable('jne',np.int32,('i',))
    latv=f.createVariable('lat',np.float64,('i',))
    lonv=f.createVariable('lon',np.float64,('i',))
    slatv=f.createVariable('start_lat',np.float64,('i',))
    slonv=f.createVariable('start_lon',np.float64,('i',))
    sdv=f.createVariable('start_d',np.float64,('i',))
    depthv=f.createVariable('depth',np.float64,('i',))
    klevv=f.createVariable('k_fixed',np.float64,('i',))
    uv=f.createVariable('uvel',np.float64,('i',))
    vv=f.createVariable('vvel',np.float64,('i',))
    dv=f.createVariable('time',np.float64,('i',))
    id_ijv=f.createVariable('id_ij',np.int32,('i',))
    id_cntv=f.createVariable('id_cnt',np.int32,('i',))

    f.file_format_major_version=1
    f.file_format_minor_version=1
    f.time_axis = 0
    iv.long_name='identification of the drifter'
    iv.units='dimensionless'
    iv.packing=0
    inev.long_name='i index'
    inev.units='none'
    inev.packing=0
    jnev.long_name='j index'
    jnev.units='none'
    jnev.packing=0
    lonv.long_name='longitude'
    lonv.units='degrees_E'
    latv.long_name='latitude'
    latv.units='degrees_N'
    slonv.long_name='start longitude'
    slonv.units='degrees_E'
    slatv.long_name='start latitude'
    slatv.units='degrees_N'
    uv.long_name= 'zonal velocity'
    uv.units='m/s'
    vv.long_name= 'meridional velocity'
    vv.units='m/s'
    sdv.long_name= 'depth of starting location'
    sdv.units='m'
    depthv.long_name= 'depth below surface'
    depthv.units='m'
    klevv.long_name= 'layer below surface'
    klevv.units='none'
    dv.units='days since 1900-01-01 00:00:00'

    id_cntv.long_name='counter component of particle id'
    id_cntv.units='dimensionless'
    id_cntv.packing=0
    id_ijv.long_name='position component of particle id'
    id_ijv.units='dimensionless'
    id_ijv.packing=0

    ivv[:]=len(ids[:])
    iv[:] = ids[:]
    inev[:]=ines[:]
    jnev[:]=jnes[:]
    lonv[:]=lons[:]
    latv[:]=lats[:]
    slonv[:]=slons[:]
    slatv[:]=slats[:]
    uv[:]=us[:]
    vv[:]=vs[:]
    sdv[:]=depths[:]
    depthv[:]=depths[:]
    klevv[:]=kfixed[:]
    dv[:]=days[:]
    f.sync()
    f.close()

lon1D = np.array([-25.,-25.,-25.,-25.,-25.,-160.,-160.,-160.,-160.])
lat1D = np.array([ 26., 26., 26., 26., 26.,   0.,   0.,   0.,   0.])
depthv= np.array([  1.,  0,   0 ,  0 ,  0 ,   0 ,   0,    0,    0 ])
klev  = np.array([ -1 ,  1 , 10 , 20 , 50 ,   1 ,  10,   20,   50 ])  # -1 if initial depth is specified, positive if keeping k level fixed
nparts=lon1D.size
lons_list = lon1D.tolist()
lats_list = lat1D.tolist()
depth_list = depthv.tolist()
klev_list = klev.tolist()
id_list = np.arange(1,nparts+1,1).tolist()
ines_list = [0] * nparts
jnes_list = [0] * nparts
days_list = [0] * nparts
u_list = [0] * nparts
v_list = [0] * nparts
id_list

createParticlesNC('drifters.res.nc',
                  ids=id_list,ines=ines_list ,jnes=jnes_list,
                  lons=lons_list,lats=lats_list,slons=lons_list,slats=lats_list,
                  us=u_list,vs=v_list,
                  depths=depth_list,kfixed=klev_list,days=days_list)





