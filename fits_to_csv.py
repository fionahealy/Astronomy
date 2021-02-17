import numpy
from matplotlib import pyplot
from astropy.io import fits
from astropy import units
from radio_beam import Beam

def fitsconverter(fitsfile):

	image = fits.getdata(fitsfile, ext=0)
	dim = len(image.shape)-2
	return image[(0,) * dim + (slice(None),)*2]



def plot_beam(minor,major,bpa,xpos,ypos,cellsize,colour,linesize):

	bpa = bpa+90	#argh why
	a=minor/(2*cellsize)     
	b=major/(2*cellsize)    
	t = numpy.linspace(0, 2*numpy.pi, 100)

	pyplot.plot(xpos+(a*numpy.cos(t)*numpy.cos(bpa) - b*numpy.sin(t)*numpy.sin(bpa)) ,\
	            ypos+(b*numpy.sin(t)*numpy.cos(bpa) + a*numpy.cos(t)*numpy.sin(bpa)) ,\
	            color=colour,linewidth = linesize)



def contourplot(image_array,levs,colour):

	peak = numpy.max(image_array)
	pyplot.contour(im,levels = [(x/100)*peak for x in levs],colors=colour,linewidths=0.5)



def raster(image_array,colour_map,origin,clim_lower,clim_upper,bar_label):
	pyplot.imshow(image_array,cmap=colour_map,origin=origin)
	pyplot.clim(-20, 20)
	cb = pyplot.colorbar()
	cb.set_label(label='RM (rad/m/m)',fontsize=12,fontname='monospace')



def rotate_grid(angle,imsize,unrotated_point):

	x_centre = imsize/2
	y_centre = imsize/2 +1

	centre = numpy.array([x_centre,y_centre])
	unrotated_point = numpy.array(unrotated_point)

	diff = unrotated_point-centre
	T = numpy.array([[numpy.cos(angle),-numpy.sin(angle)],[numpy.sin(angle),numpy.cos(angle)]])
	T_diff = numpy.dot(T,diff)
	rotated_point = T_diff + centre

	return(rotated_point)



def box(imsize,blc,trc,rot=0,pix_from_aips=True):

	# try to vectorise
	x = [blc[0],blc[0],trc[0],trc[0],blc[0]]
	y = [blc[1],trc[1],trc[1],blc[1],blc[1]]

	x_rot = []
	y_rot = []

	if rot!=0:
		for i in range(0,5):
			point_rot = rotate_grid(30,imsize,[x[i],y[i]])
			x_rot.append(point_rot[0])
			y_rot.append(point_rot[1])

		x = x_rot
		y = y_rot


	if pix_from_aips:
		y = [imsize - t for t in y]

	pyplot.plot(x,y,color='gray',linewidth=0.8,linestyle='--')



def plot_setup(title,blc,trc,bmin,bmaj,bpa,cellsize):

	corner = [blc[0]+15,blc[1]+15]

	plot_beam(bmin,bmaj,bpa,corner[0],corner[1],cellsize,"black",0.7)

	pyplot.axes().set_aspect('equal')
	pyplot.title(title,fontsize=12,fontname='monospace')

	pyplot.xlim(blc[0],trc[0])
	pyplot.ylim(blc[1],trc[1])

	pyplot.xlabel("pix",fontsize=10, fontname='monospace')
	pyplot.ylabel("pix",fontsize=10, fontname='monospace')



def get_info(file):

	info_dict = {}

	hdul = fits.open(file)
	header = hdul[0].header
	im_beam = Beam.from_fits_header(header)

	info_dict['cellsize'] = header['CDELT2']*3600
	info_dict['imsize'] = header['NAXIS1']
	info_dict['bmaj'] = im_beam.major.value*3600
	info_dict['bmin'] = im_beam.minor.value*3600
	info_dict['bpa'] = im_beam.pa.value

	return info_dict



info_dict=get_info('1800+7828I1.FITS')

im = fitsconverter('1800+7828I1.FITS')
rm = fitsconverter('1800CLIP7')

plot_setup("RM Map of 1803+784",[220,180],[400,350],info_dict['bmin'],info_dict['bmaj'],info_dict['bpa'],info_dict['cellsize'])
contourplot(im,[-0.25,0.25,0.5,1,2,4,8,16,32,64,95],'black')
raster(rm,'cubehelix','lower',-20,20,'RM (rad/m/m)')
box(info_dict['imsize'],[239,234],[275,281])
pyplot.show()


