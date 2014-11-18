
from astropy.io import fits
def read_file(filename):
    hdulist = fits.open(filename)
    hdr = hdulist[0].header
    f=open('osirisKeywordsDescription.html','w+')
    for i in hdr.cards:
        kw_name=i[0]
        kw_val = i[1]
        kw_comment=i[2]
        kw_type_raw = str(type(kw_val)).split(' ')[1]
        kw_type = kw_type_raw.replace('"','').replace("'","").replace(">","")
        f.write('<b>'+str(kw_name)+'      (native-keyword)</b><p>'+str(kw_comment)+'<p>Type='+str(kw_type)+'<p>'+'\n')
    f.close()
