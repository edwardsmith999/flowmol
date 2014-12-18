# coding: utf-8
import mayavi.mlab as mlab
import numpy as np
import sys

sys.path.insert(0,'../../../')
import postproclib as ppl

f = mlab.figure(bgcolor=(1.,1.,1.),fgcolor=(0.,0.,0.), size=(1024, 768))
#fdirs = ['/media/My Passport/Work/MD_turbulence/COMBINED_iter0_to_5600000/bins64x256x64/','/media/My Passport/Work/MD_turbulence/COMBINED_iter0_to_5600000/bins84x198x50/']
fdirs = ['/media/My Passport/Work/MD_turbulence/COMBINED_iter0_to_5600000/bins84x198x50/']
startrecs = [2805]
endrecs   = [3500]
count=548; skip = 5

#Vector component of interest
comp = 0
for i, fdir in enumerate(fdirs):
    startrec = startrecs[i] #60 #1006+93*10
    endrec = endrecs[i] #1000
    recrange = endrec - startrec + 1

    #Get initial field
    vobj = ppl.MD_vField(fdir)
    uvw = vobj.read(startrec=startrec,endrec=startrec,quit_on_error=True)
    [x,y,z] = vobj.grid
    [X,Y,Z] = np.meshgrid(x,y,z,indexing='ij')


    xh=int(np.round(uvw.shape[1]/2.))
    zh=xh
    xcell=x.max()/x.shape[0]
    ycell=y.max()/y.shape[0]
    zcell=z.max()/z.shape[0]

    avi_u = mlab.pipeline.scalar_field(X,Y,Z,uvw[:,:,:,0,comp])
    avi_uxfront = mlab.pipeline.scalar_field(X[:,:xh,:],
                                             Y[:,:xh,:],
                                             Z[:,:xh,:],
                                           uvw[:,:xh,:,0,comp])
    avi_uzfront = mlab.pipeline.scalar_field(X[:,:zh,:],
                                             Y[:,:zh,:],
                                             Z[:,:zh,:],
                                           uvw[:,:zh,:,0,comp])

    xskip = 2. ; yskip = 16.; zskip = 2.
    avi_uvw = mlab.pipeline.vector_field(X[::xskip,::yskip,::zskip], 
                                         Y[::xskip,::yskip,::zskip], 
                                         Z[::xskip,::yskip,::zskip], 
                                         uvw[::xskip,::yskip,::zskip,0,0],
                                         uvw[::xskip,::yskip,::zskip,0,1],
                                         uvw[::xskip,::yskip,::zskip,0,2])

    avi_uvwxfront = mlab.pipeline.vector_field(X[::xskip,:xh:yskip,::zskip], 
                                               Y[::xskip,:xh:yskip,::zskip], 
                                               Z[::xskip,:xh:yskip,::zskip], 
                                #np.zeros(uvw[::xskip,:xh:yskip,::zskip,0,0].shape[0:3]),
                                         uvw[::xskip,:xh:yskip,::zskip,0,0],
                                         uvw[::xskip,:xh:yskip,::zskip,0,1],
                                         uvw[::xskip,:xh:yskip,::zskip,0,2])

    avi_uvwzfront = mlab.pipeline.vector_field(X[::xskip,:zh:yskip,::zskip], 
                                               Y[::xskip,:zh:yskip,::zskip], 
                                               Z[::xskip,:zh:yskip,::zskip], 
                                         uvw[::xskip,:zh:yskip,::zskip,0,0],
                                         uvw[::xskip,:zh:yskip,::zskip,0,1],
                                         uvw[::xskip,:zh:yskip,::zskip,0,2])
                                #np.zeros(uvw[::xskip,:zh:yskip,::zskip,0,2].shape[0:3]))

    #Set camera to angle
    def camera_reset(f=f):
        mlab.view(0,180)
        f.scene.camera.azimuth(-55)
        f.scene.camera.elevation(20)
        mlab.move(0., 0., -200.)
        #f.scene.camera(0., 0., 100.)
        f.scene.camera.orthogonalize_view_up()
        f.scene.render()


    def get_removed_normal_vec(plot_vector,po):
        if po == 'x_axes':
            plot_vector_nonormal = mlab.pipeline.vector_field(
                                      plot_vector.mlab_source.x, 
                                      plot_vector.mlab_source.y, 
                                      plot_vector.mlab_source.z, 
                                      np.zeros(plot_vector.mlab_source.vectors[:,:,:,0].shape),
                                      plot_vector.mlab_source.vectors[:,:,:,1],
                                      plot_vector.mlab_source.vectors[:,:,:,2]) 
        elif po == 'y_axes':
            plot_vector_nonormal = mlab.pipeline.vector_field(
                                      plot_vector.mlab_source.x, 
                                      plot_vector.mlab_source.y, 
                                      plot_vector.mlab_source.z, 
                                      plot_vector.mlab_source.vectors[:,:,:,0],
                                      np.zeros(plot_vector.mlab_source.vectors[:,:,:,1].shape),
                                      plot_vector.mlab_source.vectors[:,:,:,2]) 
        elif po == 'z_axes':
            plot_vector_nonormal = mlab.pipeline.vector_field(
                                      plot_vector.mlab_source.x, 
                                      plot_vector.mlab_source.y, 
                                      plot_vector.mlab_source.z, 
                                      plot_vector.mlab_source.vectors[:,:,:,0],
                                      plot_vector.mlab_source.vectors[:,:,:,1],
                                      np.zeros(plot_vector.mlab_source.vectors[:,:,:,2].shape)) 
        else:
            #po is wrong, but pass and let mayavi raise the error
            plot_vector_nonormal = plot_vector
            pass

        return plot_vector_nonormal

    def slice_plane(plot_scalar=avi_u,
                    extent=[0,x.max(),
                            0,y.max(),
                            0,z.max()],
                    po='x_axes',
                    opacity = 1.,
                    xorigin = x.max()/2.,
                    yorigin = y.max()/2.,
                    zorigin = z.max()/2.):

        cp = mlab.pipeline.scalar_cut_plane( plot_scalar, 
                                             figure=f,
                                             plane_orientation=po,
                                             extent=extent,
                                             colormap='RdYlBu',
                                             opacity=1.,
                                             view_controls=False,
                                             vmin = -1.1,
                                             vmax =  1.1)

        cp.actor.property.opacity = opacity
        cp.implicit_plane.widget.origin = [xorigin, yorigin, zorigin]

        return cp

    def vector_plane(plot_vector=avi_uvw,
                    po='y_axes',
                    xorigin = x.max()/2.,
                    yorigin = y.max()/2.,
                    zorigin = z.max()/2.,
                    remove_normal=False):

    #    if remove_normal:
    #        plot_vector = get_removed_normal_vec(plot_vector,po)
    #    else:
    #        plot_vector = plot_vector

        cp = mlab.pipeline.vector_cut_plane( plot_vector, 
                                             figure=f,
                                             plane_orientation=po,
                                             view_controls=False,
                                             mode='arrow')

        cp.glyph.color_mode = 'no_coloring'
        cp.glyph.glyph.scale_mode = 'data_scaling_off'
        cp.actor.property.line_width = 4.
        cp.glyph.glyph_source.glyph_position = 'tail'
        cp.glyph.glyph.scale_factor = 50.0
        cp.glyph.glyph_source.glyph_source.tip_length = 0.2
        cp.glyph.glyph_source.glyph_source.tip_radius = 0.08
        cp.glyph.glyph_source.glyph_source.shaft_radius = 0.02
        cp.glyph.glyph_source.glyph_source.tip_resolution = 20
        cp.actor.property.opacity = 0.9

        cp.implicit_plane.widget.origin = [xorigin, yorigin, zorigin]

        return cp

    f.scene.disable_render = True

    #We seem to need a shift as no plot exactly on edges
    vshift = 14.; shift = vshift*2.
    xshift = shift*3.
    plots = []; vplots=[]

    #Front x
    plots.append(slice_plane(po='x_axes',
                             plot_scalar = avi_uxfront,
                             extent=[0,x.max(),0,y.max()/2.,shift,z.max()-vshift], 
                             xorigin=x.max()-xshift,
                             opacity = 1.))
    vplots.append(vector_plane(po='x_axes',
                               plot_vector = avi_uvwxfront,
                               xorigin=x.max()-xshift))
    #Back x
    plots.append(slice_plane(po='x_axes',
                             extent=[0,x.max(),0.,y.max(),shift,z.max()], 
                             xorigin=vshift))
    vplots.append(vector_plane(po='x_axes',
                               plot_vector = avi_uvw,
                               xorigin=shift*2.,
                               remove_normal=True))

    #Centre y
    plots.append(slice_plane(po='y_axes',
                             extent=[0.,x.max()-xshift,0,y.max(),shift-2.,z.max()]))
    vplots.append(vector_plane(po='y_axes',
                               plot_vector = avi_uvw,
                               yorigin=y.max()/2.+shift,
                               remove_normal=True))

    #Back z
    plots.append(slice_plane(po='z_axes',
                             extent=[0.,x.max()-xshift,.0,y.max(),0,z.max()], 
                             zorigin=z.max()-vshift))
    vplots.append(vector_plane(po='z_axes',
                               plot_vector = avi_uvw,
                               zorigin=z.max()-shift*2.,
                               remove_normal=True))

    #Front z
    plots.append(slice_plane(po='z_axes',
                             plot_scalar = avi_uxfront,
                             extent=[0.,x.max()-xshift,0,y.max()/2.,0,z.max()], 
                             zorigin=shift,
                             opacity = 1.))
    vplots.append(vector_plane(po='z_axes',
                               plot_vector = avi_uvwzfront,
                               zorigin=vshift))

    mlab.colorbar()
    camera_reset()


    #Loop through all values
    for rec in range(startrec,endrec,skip):

        try:
            uvw = vobj.read(startrec=rec,endrec=rec,quit_on_error=True)
            print('Reading record ', rec)
        except:
            print('Record ', rec, ' not found, skipping')
            continue

        f.scene.disable_render = True

        #Update all contour plots in figure
        for plot in plots:
            avi_u.mlab_source.set(scalars=uvw[:,:,:,0,comp])
            avi_uxfront.mlab_source.set(scalars=uvw[:,:xh,:,0,comp])
            avi_uzfront.mlab_source.set(scalars=uvw[:,:zh,:,0,comp])

        #Update all vector plots in figure
        for plot in vplots:
            avi_uvw.mlab_source.set(vectors=uvw[::xskip,::yskip,::zskip,0,:])
            avi_uvwxfront.mlab_source.set(vectors=uvw[::xskip,:xh:yskip,::zskip,0,:])
            avi_uvwzfront.mlab_source.set(vectors=uvw[::xskip,:zh:yskip,::zskip,0,:])

    #        uvwvector = np.empty((uvw.shape[0],uvw.shape[1],uvw.shape[2],uvw.shape[4]))
    #        uvwvector[:,:,:,0] = np.zeros(uvw.shape[0:3])
    #        uvwvector[:,:,:,1] = uvw[:,:,:,0,1] 
    #        uvwvector[:,:,:,2] = uvw[:,:,:,0,2]
    #        avi_uvwxfront.mlab_source.set(vectors=uvwvector[::xskip,:xh:yskip,::zskip,:])

    #        uvwvector[:,:,:,0] = uvw[:,:,:,0,0]
    #        uvwvector[:,:,:,1] = uvw[:,:,:,0,1] 
    #        uvwvector[:,:,:,2] = np.zeros(uvw.shape[0:3])
    #        avi_uvwzfront.mlab_source.set(vectors=uvwvector[::xskip,:zh:yskip,::zskip,:])

            plot.glyph.scale_mode = 'data_scaling_off'

        f.scene.disable_render = False
        filename = 'slices{num:05d}.png'.format(num=count)
        print('Wrinting', filename)
        mlab.savefig(filename)
        count += 1


#slice_plane(plane_orientation='x_axes',
#                extent=[0,x.max(),
#                        0,y.max(),
#                        0,z.max()])
#mlab.pipeline.vector_cut_plane(avi_uvw,figure=f,plane_orientation='x_axes',
#                               color=(0., 0., 0.))
#mlab.pipeline.scalar_cut_plane(avi_u,plane_orientation='y_axes',
#                                extent=[0,x.max(),0,y.max(),0,z.max()],
#                                colormap='RdYlBu',opacity=1.)
#mlab.pipeline.vector_cut_plane(avi_uvw,plane_orientation='x_axes',
#                                extent=[0,5,0,2,0,5])
#mlab.outline()
