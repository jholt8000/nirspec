
import nirspec_wavelength_utils
import numpy as np
import unittest
import nirspec_test_data

#TODO need to define dx, dy inputs for these tests to work

### nirspec_wavelength_utils tests ###
def test_2d_lambda_fit():
    """"""
    onas= [ 0.02631579,  0.02631579,  0.02631579,  0.02631579,  0.02631579,  0.02631579,
  0.02631579,  0.02631579,  0.02702703,  0.02702703,  0.02702703,  0.02702703,
  0.02702703,  0.02702703,  0.02777778,  0.02777778,  0.02777778,  0.02777778,
  0.02777778,  0.02777778,  0.02777778,  0.02777778,  0.02777778,  0.02777778,
  0.02777778,  0.02777778,  0.02777778,  0.02777778,  0.02777778,  0.02777778,
  0.02777778,  0.02777778,  0.02857143,  0.02857143,  0.02857143,  0.02857143,
  0.02941176,  0.02941176,  0.02941176,  0.02941176,  0.03030303,  0.03030303,
  0.03030303,  0.03030303]
    msls= [ 20004.971,  20008.148,  20033.158,  20068.877,  20102.977,  20115.766,
  20174.889,  20193.02,  20499.357,  20562.953,  20564.143,  20672.668,
  20728.17,   20729.859,  21033.09,   21053.775,  21061.498,  21062.689,
  21067.729,  21096.316,  21104.275,  21107.104,  21115.652,  21155.926,
  21176.387,  21232.354,  21249.48,   21278.219,  21279.893,  21316.18,
  21319.652,  21325.914,  21636.973,  21710.932,  21802.111,  21873.348,
  22311.799,  22313.65,   22460.602,  22516.727,  22938.246,  22939.947,
  22983.535,  22987.475]
    orig_pix_x_stack= [280, 291, 379, 509, 613, 663, 860, 934, 135, 356, 360, 728, 912, 918,  29,
 157, 246, 272, 281, 311, 445, 512, 695, 751, 843, 849, 965, 976, 996,  61, 299, 592, 817,
 192, 198, 657, 829,  48,  54, 190, 202]

    p1 = twod_lambda_fit.twodfit(np.array(orig_pix_x_stack),
                                         np.array(order_number_array_stack),
                                         np.array(matched_sky_line_stack), logger=self.logger,
                                         lower_len_points=10., sigma_max=0.5)

    # ## all this belongs elsewhere ###
    newdx = np.arange(1024)
    newy = 1. / self.order_num
    newoh = np.ravel(
        p1[0] + p1[1] * newdx + p1[2] * newdx ** 2 + p1[3] * newy + p1[4] * newdx * newy + p1[5] * (
        newdx ** 2) * newy)

    pl.plot(newoh,newdx)


def test_read_oh():
    # do not need dx and, just that low_d is False
    t1 = nirspec_wavelength_utils.LineId(False, [], [])
    t1.read_OH()
    if t1.ohx = nirspec_test_data.hardcoded.ohx:
        print 'nirspec_wavelength_utils.read_OH ohx passed'
    else:
        print 'nirspec_wavelength_utils.read_OH ohx failed'

    if t1.ohy = nirspec_test_data.hardcoded.ohy:
        print 'nirspec_wavelength_utils.read_OH ohy passed'
    else:
        print 'nirspec_wavelength_utils.read_OH ohy failed'

def test_gauss_sky():
    # does not need dx, dy, or low_d, perhaps this method belongs elsewhere?
    t1 = nirspec_wavelength_utils.LineId(False, [], [])
    ohx=[100,200,220.1, 300]
    ohy=[0.1, 10, 50, 74.3]
    t1.gauss_sky(ohx, ohy, sig=0.2)
    return t1.fake_sky

def test_find_xcorr_shift(ohgs):
    t1 = nirspec_wavlelength_utils.LineId(False, dx, dy)


def test_identify(self, ohx, ohy):
    pass

def test_sanity_check(orig_pix_x, order_number_array, matched_sky_line):
    pass


### nirspec_util.nirspecBookkeeping ###

def test_close_logger(self):
       pass

def test_setup_logger(self):
        """
        creates logging utility, returns logger object
        """
        log = logging.getLogger(__name__)
        x = list(log.handlers)
        for i in x:
            log.removeHandler(i)
            i.flush()
            i.close()

        if '/' in self.fits_name:
            fits_nameforlog = self.fits_name.split('/')[-1]
        else:
            fits_nameforlog = self.fits_name
        if os.path.isfile(self.outpath + '/drp_' + fits_nameforlog.rstrip('.fits') + '.log'):
            if os.path.isfile(self.outpath + '/drp_' + fits_nameforlog.rstrip('.fits') + '_old.log'):
                os.remove(self.outpath + '/drp_' + fits_nameforlog.rstrip('.fits') + '_old.log')
            os.rename(self.outpath + '/drp_' + fits_nameforlog.rstrip('.fits') + '.log',
                      self.outpath + 'drp_' + fits_nameforlog.rstrip('.fits') + '_old.log')

        logger = logging.getLogger(__name__)
        logger.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

        if not logger.handlers:
            fh = logging.FileHandler(filename=self.outpath + 'drp_' + fits_nameforlog.rstrip('.fits') + '.log')
            fh.setLevel(logging.DEBUG)
            fh.setFormatter(formatter)
            logger.addHandler(fh)

            if self.verbose:
                ch = logging.StreamHandler()
                ch.setLevel(logging.DEBUG)
                ch.setFormatter(formatter)
                logger.addHandler(ch)
            else:
                logger.propagate = False

        return logger

### nirspec_util no class
def test_make_nirspec_final_fits_and_plots(allreduceobj, order_num, sciorder, lineobj, traceobj, flatobj, sciobj):


### nirspec_util.NirspecHeader need to init a header

def test_get_data_dict(self):
    pass

def test_get_theory_order_pos(self, order_num, a_to_mu=True):
    pass

def test_sanity_check(self, flatheader):
    pass

def test_get_actual_order_num_pos(edges, theory, threshold):
    pass


#### twod_lambda_fit.py
def test_twod_lambda_fit.twodfit(dataX, dataY, dataZ, logger, lower_len_points=10., sigma_max=0.5):
    pass

### Test repo

### reduce_order

def test_reduce_order(reduction, sciobj, flatobj):
    pass


### trace_order
def __init__(self, reductionobj, flatobj, sciobj):
    self.reductionobj = reductionobj
    self.flatobj = flatobj
    self.sciobj = sciobj
    self.trace_success = False

def trace_order(self):
    pass

def test_fudge_padding(self):
    pass

def test_smooth_spectroid(self):
    pass


def test_shift_order(self):
    pass


def test_shift_order_back(self):
    pass

def test_determine_lhs_edge_pos(self, edges, theory_lhs, data_dict):
    pass