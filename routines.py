#!/usr/bin/env python
# -*- coding: utf-8 -*-
# HBV-96 model RoutinesClass

class RoutineProcess(object):
	"""docstring for RoutineProcess"""

	""""""
	
	def _precipitation(par, intab):
	    '''
	    ==============
	    Precipitation
	    ==============

	    Precipitaiton routine of the HBV96 model.

	    If temperature is lower than par['ltt'], all the precipitation is considered as
	    snow. If the temperature is higher than par['utt'], all the precipitation is
	    considered as rainfall. In case that the temperature is between par['ltt'] and
	    par['utt'], precipitation is a linear mix of rainfall and snowfall.

	    Parameters
	    ----------
	    intab['temp'] : float
	        Measured temperature [C]
	    par['ltt'] : float
	        Lower temperature treshold [C]
	    par['utt'] : float
	        Upper temperature treshold [C]
	    prec : float 
	        Precipitation [mm]
	    par['rfcf'] : float
	        Rainfall corrector factor
	    par['sfcf'] : float
	        Snowfall corrector factor

	    Returns
	    -------
	    _rf : float
	        Rainfall [mm]
	    _sf : float
	        Snowfall [mm]
	    '''

	    if intab['temp'] <= par['ltt']:
	        _rf = 0.0
	        _sf = prec*par['sfcf']

	    elif intab['temp'] >= par['utt']:
	        _rf = prec*par['rfcf']
	        _sf = 0.0

	    else:
	        _rf = ((intab['temp']-par['ltt'])/(par['utt']-par['ltt'])) * prec * par['rfcf']
	        _sf = (1.0-((intab['temp']-par['ltt'])/(par['utt']-par['ltt']))) * prec * par['sfcf']

	    return _rf, _sf


	def _snow(par, intab):
	    '''
	    ====
	    Snow
	    ====
	    
	    Snow routine of the HBV-96 model.
	    
	    The snow pack consists of two states: Water Content (wc) and Snow Pack 
	    (sp). The first corresponds to the liquid part of the water in the snow,
	    while the latter corresponds to the solid part. If the temperature is 
	    higher than the melting point, the snow pack will melt and the solid snow
	    will become liquid. In the opposite case, the liquid part of the snow will
	    refreeze, and turn into solid. The water that cannot be stored by the solid
	    part of the snow pack will drain into the soil as part of infiltration.

	    Parameters
	    ----------
	    par['cfmax'] : float 
	        Day degree factor
	    par['tfac'] : float
	        Temperature correction factor
	    intab['temp'] : float 
	        Temperature [C]
	    par['ttm'] : float 
	        Temperature treshold for Melting [C]
	    par['cfr'] : float 
	        Refreezing factor
	    par['cwh'] : float 
	        Capacity for water holding in snow pack
	    _rf : float 
	        Rainfall [mm]
	    _sf : float 
	        Snowfall [mm]
	    intab['wc'] : float 
	        Water content in previous state [mm]
	    intab['sp'] : float 
	        snow pack in previous state [mm]

	    Returns
	    -------
	    _in : float 
	        Infiltration [mm]
	    _wc_new : float 
	        Water content in posterior state [mm]
	    _intab['sp'] : float 
	        Snowpack in posterior state [mm]
	    '''

	    if intab['temp'] > par['ttm']:

	        if par['cfmax']*(intab['temp']-par['ttm']) < intab['sp']+_sf:
	            _melt = par['cfmax']*(intab['temp']-par['ttm'])
	        else:
	            _melt = intab['sp']+_sf

	        _intab['sp'] = intab['sp'] + _sf - _melt
	        _wc_int = intab['wc'] + _melt + _rf

	    else:
	        if par['cfr']*par['cfmax']*(par['ttm']-intab['temp']) < intab['wc']:
	            _refr = par['cfr']*par['cfmax']*(par['ttm'] - intab['temp'])
	        else:
	            _refr = intab['wc'] + _rf

	        _intab['sp'] = intab['sp'] + _sf + _refr
	        _wc_int = intab['wc'] - _refr + _rf

	    if _wc_int > par['cwh']*_intab['sp']:
	        _in = _wc_int-par['cwh']*_intab['sp']
	        _wc_new = par['cwh']*_intab['sp']
	    else:
	        _in = 0.0
	        _wc_new = _wc_int

	    return _in, _wc_new, _intab['sp']

	
	def _soil(par, intab):
	    '''
	    ====
	    Soil
	    ====
	    
	    Soil routine of the HBV-96 model.
	    
	    The model checks for the amount of water that can infiltrate the soil, 
	    coming from the liquid precipitation and the snow pack melting. A part of 
	    the water will be stored as soil moisture, while other will become runoff, 
	    and routed to the upper zone tank.

	    Parameters
	    ----------
	    par['fc'] : float 
	        Filed capacity
	    par['beta'] : float 
	        Shape coefficient for effective precipitation separation
	    par['etf'] : float 
	        Total potential evapotranspiration
	    intab['temp'] : float 
	        Temperature
	    intab['tm'] : float 
	        Average long term temperature
	    par['e_corr'] : float 
	        Evapotranspiration corrector factor
	    par['lp'] : float _soil 
	        wilting point
	    par['tfac'] : float 
	        Time conversion factor
	    par['c_flux'] : float 
	        Capilar flux in the root zone
	    _in : float 
	        actual infiltration
	    intab['ep'] : float 
	        actual evapotranspiration
	    intab['sm'] : float 
	        Previous soil moisture value
	    intab['uz'] : float 
	        Previous Upper zone value

	    Returns
	    -------
	    intab['sm'] : float 
	        New value of soil moisture
	    intab['uz'] : float 
	        New value of direct runoff into upper zone
	    '''

	    qdr = max(intab['sm'] + inf - par['fc'], 0)
	    _in = inf - qdr
	    _r = ((intab['sm']/par['fc'])**par['beta']) * _in
	    _ep_int = (1.0 + par['etf']*(intab['temp'] - intab['tm']))*par['e_corr']*intab['ep']
	    _ea = max(_ep_int, (intab['sm']/(par['lp']*par['fc']))*_ep_int)

	    _cf = par['c_flux']*((par['fc'] - intab['sm'])/par['fc'])
	    intab['sm'] = max(intab['sm'] + _in - _r + _cf - _ea, 0)
	    intab['uz'] = intab['uz'] + _r - _cf

	    return intab['sm'], intab['uz'], qdr


	def _response(par, intab):
	    '''
	    ========
	    Response
	    ========
	    The response routine of the HBV-96 model.
	    
	    The response routine is in charge of transforming the current values of 
	    upper and lower zone into discharge. This routine also controls the 
	    recharge of the lower zone tank (baseflow). The transformation of units 
	    also occurs in this point.
	    
	    Parameters
	    ----------
	    par['tfac'] : float
	        Number of hours in the time step
	    par['perc'] : float
	        Percolation value [mm\hr]
	    par['alpha'] : float
	        Response box parameter
	    par['k'] : float
	        Upper zone recession coefficient
	    par['k1'] : float 
	        Lower zone recession coefficient
	    par['area'] : float
	        Catchment par['area'] [Km2]
	    intab['lz'] : float 
	        Previous lower zone value [mm]
	    intab['uz'] : float 
	        Previous upper zone value before percolation [mm]
	    qdr : float
	        Direct runoff [mm]
	    
	    '''    
	    intab['lz'] = intab['lz'] + np.min(par['perc'], intab['uz'])
	    intab['uz'] = np.max(intab['uz'] - par['perc'], 0.0)

	    _q_0 = par['k']*(intab['uz']**(1.0 + par['alpha']))
	    _q_1 = par['k1']*intab['lz']

	    intab['uz'] = max(intab['uz'] - (_q_0), 0)
	    intab['lz'] = max(intab['lz'] - (_q_1), 0)

	    q_new = par['area']*(_q_0 + _q_1 + qdr)/(3.6*par['tfac'])

	    return q_new, intab['uz'], intab['lz']


	def _routing():
	    """
	    To be implemented
	    """
	    return


	def _step_run(par, intab):
	    '''
	    ========
	    Step run
	    ========
	    
	    Makes the calculation of next step of discharge and states
	    
	    Parameters
	    ----------
	    p : array_like [18]
	        Parameter vector, set up as:
	        [par['ltt'], par['utt'], par['ttm'], par['cfmax'], par['fc'], ecorr, par['etf'], par['lp'], par['k'], par['k1'], 
	        par['alpha'], par['beta'], par['cwh'], par['cfr'], par['c_flux'], par['perc'], par['rfcf'], par['sfcf']]
	    p2 : array_like [2]
	        Problem parameter vector setup as:
	        [par['tfac'], par['area']]
	    v : array_like [4]
	        Input vector setup as:
	        [prec, intab['temp'], evap, llt]
	    St : array_like [5]
	        Previous model states setup as:
	        [sp, sm, uz, lz, wc]

	    Returns
	    -------
	    q_new : float
	        Discharge [m3/s]
	    St : array_like [5]
	        Posterior model states
	    '''
	    ## Parse of parameters from input vector to model
	    # par['ltt'] = p[0]
	    # par['utt'] = p[1]
	    # par['ttm'] = p[2]
	    # par['cfmax'] = p[3]
	    # par['fc'] = p[4]
	    # par['e_corr'] = p[5]
	    # par['etf'] = p[6]
	    # par['lp'] = p[7]
	    # par['k'] = p[8]
	    # par['k1'] = p[9]
	    # par['alpha'] = p[10]
	    # par['beta'] = p[11]
	    # par['cwh'] = p[12]
	    # par['cfr'] = p[13]
	    # par['c_flux'] = p[14]
	    # par['perc'] = p[15]
	    # par['rfcf'] = p[16]
	    # par['sfcf'] = p[17]

	    # ## Non optimisable parameters
	    # par['tfac'] = p2[0]
	    # par['area'] = p2[1]

	    # ## Parse of Inputs
	    # intab['avg_prec'] = v[0] # Precipitation [mm]
	    # intab['temp'] = v[1] # Temperature [C]
	    # intab['ep'] = v[2] # Long terms (monthly) Evapotranspiration [mm]
	    # intab['tm'] = v[3] #Long term (monthly) average temperature [C]

	    # ## Parse of states
	    # intab['sp'] = St[0]
	    # intab['sm'] = St[1]
	    # intab['uz'] = St[2]
	    # intab['lz'] = St[3]
	    # intab['wc'] = St[4]

	    rf, sf = _precipitation(par, intab)
	    
	    inf, intab['wc'], intab['sp'] = _snow(par, intab)
	    
	    intab['sm'], intab['uz'], qdr = _soil(par, intab)
	    
	    q_new, intab['uz'], intab['lz'] = _response(par, intab)

	    return intab


	def simulate(intab['avg_prec'], intab['temp'], et, par, p2, init_st=DEF_ST, ll_temp=None, 
	             q_0=DEF_q0):
	    '''
	    ========
	    Simulate
	    ========

	    Run the HBV model for the number of steps (n) in precipitation. The
	    resluts are (n+1) simulation of discharge as the model calculates step n+1
		
	    
	    Parameters
	    ----------
	    intab['avg_prec'] : array_like [n]
	        Average precipitation [mm/h]
	    intab['temp'] : array_like [n]
	        Average temperature [C]
	    et : array_like [n]
	        Potential Evapotranspiration [mm/h]
	    par : array_like [18]
	        Parameter vector, set up as:
	        [par['ltt'], par['utt'], par['ttm'], par['cfmax'], par['fc'], ecorr, par['etf'], par['lp'], par['k'], par['k1'], 
	        par['alpha'], par['beta'], par['cwh'], par['cfr'], par['c_flux'], par['perc'], par['rfcf'], par['sfcf']]
	    p2 : array_like [2]
	        Problem parameter vector setup as:
	        [par['tfac'], par['area']]
	    init_st : array_like [5], optional
	        Initial model states, [sp, sm, uz, lz, wc]. If unspecified, 
	        [0.0, 30.0, 30.0, 30.0, 0.0] mm
	    ll_temp : array_like [n], optional
	        Long term average temptearature. If unspecified, calculated from intab['temp'].
	    q_0 : float, optional
	        Initial discharge value. If unspecified set to 10.0
	    

	    Returns
	    -------
	    q_sim : array_like [n]
	        Discharge for the n time steps of the precipitation vector [m3/s]
	    st : array_like [n, 5]
	        Model states for the complete time series [mm]
	    '''


	    st = [init_st, ]

	    if ll_temp is None:
	        ll_temp = [np.mean(intab['temp']), ] * len(intab['avg_prec'])

	    q_sim = [q_0, ]

	    for i in range(len(intab['avg_prec'])):
	        v = [intab['avg_prec'][i], intab['temp'][i], et[i], ll_temp[i]]
	        q_out, st_out = _step_run(par, p2, v, st[i])
	        q_sim.append(q_out)
	        st.append(st_out)

	    return q_sim, st


	def _nse(q_rec, q_sim):
	    '''
	    ===
	    NSE
	    ===
	    
	    Nash-Sutcliffe efficiency. Metric for the estimation of performance of the 
	    hydrological model
	    
	    Parameters
	    ----------
	    q_rec : array_like [n]
	        Measured discharge [m3/s]
	    q_sim : array_like [n] 
	        Simulated discharge [m3/s]
	        
	    Returns
	    -------
	    f : float
	        NSE value
	    '''
	    a = np.square(np.subtract(q_rec, q_sim))
	    b = np.square(np.subtract(q_rec, np.nanmean(q_rec)))
	    if a.any < 0.0:
	        return(np.nan)
	    f = 1.0 - (np.nansum(a)/np.nansum(b))
	    return f


	def _rmse(q_rec,q_sim):
	    '''
	    ====
	    RMSE
	    ====
	    
	    Root Mean Squared Error. Metric for the estimation of performance of the 
	    hydrological model.
	    
	    Parameters
	    ----------
	    q_rec : array_like [n]
	        Measured discharge [m3/s]
	    q_sim : array_like [n] 
	        Simulated discharge [m3/s]
	        
	    Returns
	    -------
	    f : float
	        RMSE value
	    '''
	    erro = np.square(np.subtract(q_rec,q_sim))
	    if erro.any < 0:
	        return(np.nan)
	    f = np.sqrt(np.nanmean(erro))
	    return f

	def calibrate(flow, intab['avg_prec'], intab['temp'], et, p2, init_st=None, ll_temp=None,
	              x_0=None, x_lb=P_LB, x_ub=P_UB, obj_fun=_rmse, wu=10,
	              verbose=False, tol=0.001, minimise=True, fun_nam='RMSE'):
	    '''
	    =========
	    Calibrate
	    =========

	    Run the calibration of the HBV-96. The calibration is used to estimate the
	    optimal set of parameters that minimises the difference between 
	    observations and modelled discharge.
	    
	    Parameters
	    ----------
	    
	    flow : array_like [n]
	        Measurements of discharge [m3/s]
	    intab['avg_prec'] : array_like [n]
	        Average precipitation [mm/h]
	    intab['temp'] : array_like [n]
	        Average temperature [C]
	    et : array_like [n]
	        Potential Evapotranspiration [mm/h] 
	    p2 : array_like [2]
	        Problem parameter vector setup as:
	        [par['tfac'], par['area']]
	    init_st : array_like [5], optional
	        Initial model states, [sp, sm, uz, lz, wc]. If unspecified, 
	        [0.0, 30.0, 30.0, 30.0, 0.0] mm
	    ll_temp : array_like [n], optional
	        Long term average temptearature. If unspecified, calculated from intab['temp'].
	    x_0 : array_like [18], optional
	        First guess of the parameter vector. If unspecified, a random value
	        will be sampled between the boundaries of the 
	    x_lb : array_like [18], optional
	        Lower boundary of the parameter vector. If unspecified, a random value
	        will be sampled between the boundaries of the 
	    x_ub : array_like [18], optional
	        First guess of the parameter vector. If unspecified, a random value
	        will be sampled between the boundaries of the
	    obj_fun : function, optional
	        Function that takes 2 parameters, recorded and simulated discharge. If
	        unspecified, RMSE is used.
	    wu : int, optional
	        Warming up period. This accounts for the number of steps that the model
	        is run before calculating the performance function.
	    verbose : bool, optional
	        If True, displays the result of each model evaluation when performing
	        the calibration of the hydrological model.
	    tol : float, optional
	        Determines the tolerance of the solutions in the optimisaiton process.
	    minimise : bool, optional
	        If True, the optimisation corresponds to the minimisation of the 
	        objective function. If False, the optimial of the objective function is
	        maximised.
	    fun_nam : str, optional
	        Name of the objective function used in calibration. If unspecified, is
	        'RMSE'
	    
	    Returns
	    -------
	    params : array_like [18]
	        Optimal parameter set
	    
	    performance : float
	        Optimal value of the objective function
	    '''

	    def _cal_fun(par):
	        q_sim = simulate(intab['avg_prec'][:-1], intab['temp'], et, par, p2, init_st=None,
	                         ll_temp=None, q_0=10.0)[0]
	        if minimise:
	            perf = obj_fun(flow[wu:], q_sim[wu:])
	        else:
	            perf = -obj_fun(flow[wu:], q_sim[wu:])

	        if verbose:
	            print('{0}: {1}'.format(fun_nam, perf))
	        return perf

	    # Boundaries
	    x_b = zip(x_lb, x_ub)

	    # initial guess
	    if x_0 is None:
	        # Randomly generated
	        x_0 = np.random.uniform(x_lb, x_ub)

	    # Model optimisation
	    par_cal = opt.minimize(_cal_fun, x_0, method='L-BFGS-B', bounds=x_b,
	                           tol=tol)
	    params = par_cal.x
	    performance = par_cal.fun
	    return params, performance

