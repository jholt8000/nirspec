
import wavelength_utils
reload(wavelength_utils)
import numpy as np
import unittest
import array_manipulate
reload(array_manipulate)
import keck_fits
import nirspec_util
reload(nirspec_util)
import reduce_order
reload(reduce_order)
import nirspecOO
reload(nirspecOO)
import twod_lambda_fit
reload(twod_lambda_fit)
from astropy.io import fits

class Unit_tests(unittest.TestCase):

    def __init__(self, sci_name='NS.20120929.34991.fits', flat_name='NS.20120929.09078.fits'):
        '''
         make input objects to unit tests
        :param sci_name:
        :param flat_name:
        :return:

        '''
        self.nb = nirspec_util.NirspecBookkeeping(sci_name, outpath='.')
        self.sci, self.sciheader, self.sciname = keck_fits.Handle_fits.get_array_and_header(sci_name)

        self.flat, self.flatheader, self.flatname = keck_fits.Handle_fits.get_array_and_header(flat_name)

        # instantiate sciObj and flatObj objects
        self.sciObj = array_manipulate.SciArray(self.sci)
        self.flatObj = array_manipulate.FlatArray(self.flat)

        self.reduction = nirspecOO.Main(sci_name, flat_name)

        self.reduced_order_object = reduce_order.Reduce_order(self.reduction, self.sciObj, self.flatObj)
        self.reduced_order_object.reduce_order()

        #this needs to be writing this array somewhere before the tests are run, otherwise it is testing itself
        self.setup_initial_cleaned_data = fits.getdata('cleaned2.fits')

        self.onas= np.array([ 0.02631579,  0.02631579,  0.02631579,  0.02631579,  0.02631579,  0.02631579,
          0.02631579,  0.02631579,  0.02702703,  0.02702703,  0.02702703,  0.02702703,
          0.02702703,  0.02702703,  0.02777778,  0.02777778,  0.02777778,  0.02777778,
          0.02777778,  0.02777778,  0.02777778,  0.02777778,  0.02777778,  0.02777778,
          0.02777778,  0.02777778,  0.02777778,  0.02777778,  0.02777778,  0.02777778,
          0.02777778,  0.02777778,  0.02857143,  0.02857143,  0.02857143,  0.02857143,
          0.02941176,  0.02941176,  0.02941176,  0.02941176,  0.03030303,  0.03030303,
          0.03030303,  0.03030303])
        self.msls= np.array([ 20004.971,  20008.148,  20033.158,  20068.877,  20102.977,  20115.766,
          20174.889,  20193.02,  20499.357,  20562.953,  20564.143,  20672.668,
          20728.17,   20729.859,  21033.09,   21053.775,  21061.498,  21062.689,
          21067.729,  21096.316,  21104.275,  21107.104,  21115.652,  21155.926,
          21176.387,  21232.354,  21249.48,   21278.219,  21279.893,  21316.18,
          21319.652,  21325.914,  21636.973,  21710.932,  21802.111,  21873.348,
          22311.799,  22313.65,   22460.602,  22516.727,  22938.246,  22939.947,
          22983.535,  22987.475])
        self.orig_pix_x_stack= np.array([280, 291, 379, 509, 613, 663, 860, 934, 135, 356, 360, 728, 912, 918,  29,
         157, 246, 272, 281, 311, 445, 512, 695, 751, 843, 849, 965, 976, 996,  61, 299, 592, 817,
         192, 198, 657, 829,  48,  54, 190, 202, 300, 400, 500])

        self.p1 = np.array([217.748988342,   1.48988109069,  -0.000704920991226,   749152.726451, -46.9689970409,   .0278160250772])

        self.newoh=np.array([])
        self.ohx,self.ohy=([],[])
        self.fake_sky=[]
        self.id_tuple=[]

    def run_all(self):

        #np.testing.assert_array_almost_equal(self.test_cosmic(), self.setup_initial_cleaned_data)

        p1_result = self.test_twod_lambda_fit_twodfit()
        np.testing.assert_array_almost_equal(p1_result, self.p1)

        fake_sky = self.test_gauss_sky()
        return fake_sky
        #self.assertEqual(self.test_twod_lambda_fit_applySolution(), self.newoh)
        #self.assertEqual(self.test_read_oh(), (self.ohx,self.ohy))
        #self.assertEqual(self.test_gauss_sky(),self.fake_sky)
        #self.assertEqual(self.test_find_xcorr_shift(), self.lambda_shift)
        #self.assertEqual(self.test_identify(), self.id_tuple)

    def make_cosmic_data(self):
        '''only run once when making the tests'''
        cc = self.test_cosmic()
        hdu = fits.PrimaryHDU(cc)
        hdu.writeto('cleaned.fits', clobber=True)
        cc2 = fits.getdata('cleaned.fits')
        np.testing.assert_array_almost_equal(cc, cc2 )

    def test_cosmic(self):
            self.sciObj.cosmic()
            return self.sciObj.data

    def make_twod_lamdba_fit_data(self):
        ''' only run once when making the tests'''
        p1 = self.test_twod_lambda_fit_twodfit()
        print p1

    def test_twod_lambda_fit_twodfit(self):
        test_p1, newoh, dataZZ = twod_lambda_fit.twodfit(self.orig_pix_x_stack, self.onas, self.msls, logger=self.reduction.logger,
                                             lower_len_points=10., sigma_max=0.5)
        print 'test_p1=',test_p1
        return test_p1

    def test_twod_lambda_fit_applySolution(self):
        newoh = twod_lambda_fit.applySolution(order_object, p1)
        return newoh

    def test_read_oh():
        # do not need dx and, just that low_d is False
        t1 = wavelength_utils.LineId(False, [], [])
        ohx, ohy = t1.read_OH()
        return ohx,ohy

    def test_gauss_sky(self):
        # does not need dx, dy, or low_d, perhaps this method belongs elsewhere?
        t1 = wavelength_utils.LineId(self.dx, self.sky, self.low_disp, self.logger)
        fake_sky = t1.gauss_sky(self.ohx, self.ohy, 0.2)
        return fake_sky

    def test_find_xcorr_shift(self):
        t1 = wavlelength_utils.LineId(self.dx, self.sky, self.low_disp, self.logger)
        lambda_shift = t1.find_xcorr_shift(self.fake_sky)
        return lambda_shift

    def test_identify(self):
           id_tuple = self.identify(self.ohx, self.ohy)
           return id_tuple
           #      self.matchesdx, self.matchesohx, self.matchesohy, self.bigohx, self.bigohy, \
           #     self.identify_status, self.matchesidx = id_tuple


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
#def test_make_nirspec_final_fits_and_plots(allreduceobj, order_num, sciorder, lineobj, traceobj, flatobj, sciobj):


### nirspec_util.NirspecHeader need to init a header
#
# def test_get_data_dict(self):
#     pass
#
# def test_get_theory_order_pos(self, order_num, a_to_mu=True):
#     pass
#
# def test_sanity_check(self, flatheader):
#     pass
#
# def test_get_actual_order_num_pos(edges, theory, threshold):
#     pass
#
# ### reduce_order
#
# def test_reduce_order(reduction, sciobj, flatobj):
#     pass
#
#
# ### trace_order
# def __init__(self, reductionobj, flatobj, sciobj):
#     self.reductionobj = reductionobj
#     self.flatobj = flatobj
#     self.sciobj = sciobj
#     self.trace_success = False
#
# def trace_order(self):
#     pass
#
# def test_fudge_padding(self):
#     pass
#
# def test_smooth_spectroid(self):
#     pass
#
#
# def test_shift_order(self):
#     pass
#
#
# def test_shift_order_back(self):
#     pass
#
# def test_determine_lhs_edge_pos(self, edges, theory_lhs, data_dict):
#     pass