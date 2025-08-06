import glob
import json
import matplotlib.pyplot as plt
import seaborn as sns; sns.set_context('talk')
from pycorn import pc_res3, pc_uni6
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import linregress
from scipy.signal import find_peaks
import os


def gaussian(x, a1, b1, c1):
    '''
    Single Gaussian function (1 kernel).
    '''
    
    return a1*np.exp((-(x-b1)**2)/(2*c1**2)) 

def gaussians6(x, 
               a1, a2, a3, a4, a5, a6, 
               b1, b2, b3, b4, b5, b6, 
               c1, c2, c3, c4, c5, c6):
    '''
    Explicit mixed Gaussian function (5 kernels).
    '''
    
    return ( a1*np.exp((-(x-b1)**2)/(2*c1**2)) + 
           a2*np.exp((-(x-b2)**2)/(2*c2**2)) + 
           a3*np.exp((-(x-b3)**2)/(2*c3**2)) + 
           a4*np.exp((-(x-b4)**2)/(2*c4**2)) + 
           a5*np.exp((-(x-b5)**2)/(2*c5**2))  + 
           a6*np.exp((-(x-b6)**2)/(2*c6**2))
           )        

def gaussians5(x, a1, a2, a3, a4, a5, b1, b2, b3, b4, b5, c1, c2, c3, c4, c5):
    '''
    Explicit mixed Gaussian function (5 kernels).
    '''
    
    return a1*np.exp((-(x-b1)**2)/(2*c1**2)) + a2*np.exp((-(x-b2)**2)/(2*c2**2)) + a3*np.exp((-(x-b3)**2)/(2*c3**2)) + a4*np.exp((-(x-b4)**2)/(2*c4**2)) + a5*np.exp((-(x-b5)**2)/(2*c5**2)) 



# calibrations
def S75_5_150_calibration_agilent(path_dextran, path_calib, flowrate, comment="" ):
    '''
    calibrate with 5 peaks
    '''
    
    savename = 'calibrations/' + "HPLC_" + path_dextran.split('/')[-1].replace(".CSV", "").replace(" ", "_") + "_" + comment
    

    hmw_std = {
        'Conalbumin':75000,
        'Ovalbumin':43000,
        'Carbonic anhydrase':29000,
        'Ribonuclease A':13700,
        'Aprotinin':6500
    }

    dextran = np.loadtxt(path_dextran, delimiter=",")
    dextran[:,0] = dextran[:,0] * flowrate

    v_dex, uv_dex = dextran[:,0] , dextran[:,1]
    void = v_dex[np.argmax(uv_dex)]

    print ( f"VOID: {void:1.3f} mL")

    proteins = np.loadtxt(path_calib, delimiter=",")
    proteins[:,0] = proteins[:,0] * flowrate
    v_prot, uv_prot = proteins[:,0], proteins[:,1]


    # Get peaks by fitting to Gaussian functions
    peak_pos_guess = np.array([ #void*0.95, #void peak 
                               1.22892861 , 1.33212932, 1.54220434, 1.78828965, 2.12824897,
                               #1.25, 1.35, 1.6, 1.80, 2.15
                              ])
    bounds = np.array( 
                       [ [0, 1100], ]   *len(peak_pos_guess) +     #scale
                       [ [1.0,2.6], ]   *len(peak_pos_guess) +     # mean
                       [ [1e-5, 1.0e-1], ] *len(peak_pos_guess) ,     # std
                         ).T
    p0 = np.hstack([np.array([50]*len(peak_pos_guess) ), np.array(peak_pos_guess) , np.array([0.01]*len(peak_pos_guess) )])
    print(p0.shape)
    print(bounds.shape)
    
    # popt= p0
    popt, pcov = curve_fit(gaussians5, 
                           v_prot, 
                           uv_prot, 
                           p0=p0,
                           bounds=bounds ,
                           # maxfev=1e6
                          )
 #    popt = [1.26517452e+02 ,5.04075172e+02, 5.99687981e+02, 2.83442026e+02,
 # 4.16019723e+02, 3.62879603e+02, 1.11775545e+00 ,1.22566036e+00,
 # 1.33943519e+00, 1.53629833e+00 ,1.78834227e+00 ,2.09254408e+00,
 # 4.08184419e-02, 4.26238299e-02 ,6.22584807e-02, 5.47710286e-02,
 # 5.06132447e-02, 5.60566319e-02,]
    
    print(popt)
    peaks = popt[len(peak_pos_guess)+1:len(peak_pos_guess)+len(peak_pos_guess)]
    print(peaks)
    

    # Plot
    fig, ax = plt.subplots(ncols=2, figsize=(12,5))
    ax[0].plot(v_dex, uv_dex, label='Dextran')
    ax[0].plot(v_prot, uv_prot, 'k-', label='Proteins')
    ax[0].plot(v_prot, gaussians5(v_prot, *popt), 'r-', linewidth=1, label='Fit')
    
    for i in range(len(peak_pos_guess)):
        color='r'
        if i == 0:
            color='gray'
        ax[0].fill_between(v_prot, gaussian(v_prot, *popt[i::len(peak_pos_guess)]), color=color, alpha=0.2)

        
    
    ax[0].set(xlim=[0.8,2.5])
    ax[0].legend(loc='center left', bbox_to_anchor=(1, 0.5) )#loc='best')
    ax[0].set(xlabel='Retention volume / mL', ylabel='A280 / a.u.',)# title=chroma.split('/')[-2])

    
    # Calibration curve
    Kav = (peaks - void) / (3-void)
    log10mw = np.log10(list(hmw_std.values()))
    slope, intercept, rvalue, pvalue, stderr = linregress(log10mw, Kav)
    xs = np.linspace(log10mw.min(), log10mw.max(), 100)
    
    for peak in peaks:
        ax[0].axvline(x=peak, color='gray', zorder=-10, alpha=0.5)
        ax[0].text(peak, 1.2*np.max(uv_prot), f"{peak:2.2f} mL", ha="left", va="center", rotation=60, fontsize=10)
    
    ax[0].axvline(x=void, color='k', zorder=-10, alpha=0.8)
    ax[0].text(void, 1.2*np.max(uv_prot), f"Void\n{void:2.2f} mL", ha="right", va="center", rotation=60, fontsize=10)
    
    ax[0].set_ylim( [-1.0, 1.2 * np.max(uv_prot) ] )


    ax[1].scatter(log10mw, Kav, color='white', edgecolor='k')
    ax[1].plot(xs, intercept + slope*xs, 'k-', linewidth=1, label=f'$R^2$ = {rvalue**2:.3f}')
    ax[1].legend(loc='best' ) #loc='best',
    ax[1].set(xlabel='log10(MW)', ylabel='Kav',)# title=chroma.split('/')[-2])
    plt.tight_layout()
    plt.savefig(savename + '.png', dpi=300)
    plt.show()
    
    # Save calibration data.
    calibration_data = {
        'log10mw':log10mw.tolist(),
        'Kav':Kav.tolist(),
        'Vc':3,
        'Vo':void,
        'intercept':intercept,
        'slope':slope
    }
    
    print("saving to: " , os.path.abspath ( savename + f'_excl_void_agilent.json' )  )
    with open(savename + f'_excl_void_agilent.json' , 'w') as f:
        json.dump(calibration_data, f)

def S75_5_150_calibration_agilent_fitvoid(path_dextran, path_calib, flowrate, comment="" ):
    '''
    calibrate with 6 peaks, to exclude void
    '''
    
    savename = 'calibrations/' + "HPLC_" + path_dextran.split('/')[-1].replace(".CSV", "").replace(" ", "_") + "_" + comment
    

    hmw_std = {
        'Conalbumin':75000,
        'Ovalbumin':43000,
        'Carbonic anhydrase':29000,
        'Ribonuclease A':13700,
        'Aprotinin':6500
    }

    dextran = np.loadtxt(path_dextran, delimiter=",")
    dextran[:,0] = dextran[:,0] * flowrate

    v_dex, uv_dex = dextran[:,0] , dextran[:,1]
    void = v_dex[np.argmax(uv_dex)]

    print ( f"VOID: {void:1.3f} mL")

    proteins = np.loadtxt(path_calib, delimiter=",")
    proteins[:,0] = proteins[:,0] * flowrate
    v_prot, uv_prot = proteins[:,0], proteins[:,1]


    # Get peaks by fitting to Gaussian functions
    peak_pos_guess = np.array([void*0.95, #void peak 
                               1.22892861 , 1.33212932, 1.54220434, 1.78828965, 2.12824897,
                               #1.25, 1.35, 1.6, 1.80, 2.15
                              ])
    bounds = np.array( 
                       [ [0, 1100], ]   *len(peak_pos_guess) +     #scale
                       [ [1.0,2.6], ]   *len(peak_pos_guess) +     # mean
                       [ [1e-5, 1.0e-1], ] *len(peak_pos_guess) ,     # std
                         ).T
    p0 = np.hstack([np.array([50]*len(peak_pos_guess) ), np.array(peak_pos_guess) , np.array([0.01]*len(peak_pos_guess) )])
    print(p0.shape)
    print(bounds.shape)
    
    # popt= p0
    popt, pcov = curve_fit(gaussians6, 
                           v_prot, 
                           uv_prot, 
                           p0=p0,
                           bounds=bounds ,
                           # maxfev=1e6
                          )
 #    popt = [1.26517452e+02 ,5.04075172e+02, 5.99687981e+02, 2.83442026e+02,
 # 4.16019723e+02, 3.62879603e+02, 1.11775545e+00 ,1.22566036e+00,
 # 1.33943519e+00, 1.53629833e+00 ,1.78834227e+00 ,2.09254408e+00,
 # 4.08184419e-02, 4.26238299e-02 ,6.22584807e-02, 5.47710286e-02,
 # 5.06132447e-02, 5.60566319e-02,]
    
    print(popt)
    peaks = popt[len(peak_pos_guess)+1:len(peak_pos_guess)+len(peak_pos_guess)]
    print(peaks)
    

    # Plot
    fig, ax = plt.subplots(ncols=2, figsize=(12,5))
    ax[0].plot(v_dex, uv_dex, label='Dextran')
    ax[0].plot(v_prot, uv_prot, 'k-', label='Proteins')
    ax[0].plot(v_prot, gaussians6(v_prot, *popt), 'r-', linewidth=1, label='Fit')
    
    for i in range(len(peak_pos_guess)):
        color='r'
        if i == 0:
            color='gray'
        ax[0].fill_between(v_prot, gaussian(v_prot, *popt[i::len(peak_pos_guess)]), color=color, alpha=0.2)

        
    
    ax[0].set(xlim=[0.8,2.5])
    ax[0].legend(loc='center left', bbox_to_anchor=(1, 0.5) )#loc='best')
    ax[0].set(xlabel='Retention volume / mL', ylabel='A280 / a.u.',)# title=chroma.split('/')[-2])

    
    # Calibration curve
    Kav = (peaks - void) / (3-void)
    log10mw = np.log10(list(hmw_std.values()))
    slope, intercept, rvalue, pvalue, stderr = linregress(log10mw, Kav)
    xs = np.linspace(log10mw.min(), log10mw.max(), 100)
    
    for peak in peaks:
        ax[0].axvline(x=peak, color='gray', zorder=-10, alpha=0.5)
        ax[0].text(peak, 1.2*np.max(uv_prot), f"{peak:2.2f} mL", ha="left", va="center", rotation=60, fontsize=10)
    
    ax[0].axvline(x=void, color='k', zorder=-10, alpha=0.8)
    ax[0].text(void, 1.2*np.max(uv_prot), f"Void\n{void:2.2f} mL", ha="right", va="center", rotation=60, fontsize=10)
    
    ax[0].set_ylim( [-1.0, 1.2 * np.max(uv_prot) ] )


    ax[1].scatter(log10mw, Kav, color='white', edgecolor='k')
    ax[1].plot(xs, intercept + slope*xs, 'k-', linewidth=1, label=f'$R^2$ = {rvalue**2:.3f}')
    ax[1].legend(loc='best' ) #loc='best',
    ax[1].set(xlabel='log10(MW)', ylabel='Kav',)# title=chroma.split('/')[-2])
    plt.tight_layout()
    plt.savefig(savename + '.png', dpi=300)
    plt.show()
    
    # Save calibration data.
    calibration_data = {
        'log10mw':log10mw.tolist(),
        'Kav':Kav.tolist(),
        'Vc':3,
        'Vo':void,
        'intercept':intercept,
        'slope':slope
    }
    
    print("saving to: " , os.path.abspath ( savename + f'_excl_void_agilent.json' )  )
    with open(savename + f'_excl_void_agilent.json' , 'w') as f:
        json.dump(calibration_data, f)

def S200_5_150_calibration_agilent_fitvoid(path_dextran, path_calib, flowrate, comment="" ):
    '''
    fit 6 peaks to exclude void
    '''
    
    savename = 'calibrations/' + "HPLC_" + path_dextran.split('/')[-1].replace(".CSV", "").replace(" ", "_") + "_" + comment
    
    hmw_std = {
        'Thyroglobulin':669000,
        'Ferritin':440000,
        'Aldolase':158000,
        'Conalbumin':75000,
        'Ovalbumin':43000,
    }


    dextran = np.loadtxt(path_dextran, delimiter=",")
    dextran[:,0] = dextran[:,0] * flowrate

    v_dex, uv_dex = dextran[:,0] , dextran[:,1]
    void = v_dex[np.argmax(uv_dex)]


    
    print ( f"VOID: {void:1.3f} mL")

    proteins = np.loadtxt(path_calib, delimiter=",")
    proteins[:,0] = proteins[:,0] * flowrate
    v_prot, uv_prot = proteins[:,0], proteins[:,1]

    # plt.plot(dextran[:,0], dextran[:,1])
    # plt.plot(proteins[:,0], proteins[:,1])
    # plt.show()

    # Get peaks by fitting to Gaussian functions
    peak_pos_guess = np.array([void*0.95, #void peak 
                               1.22892861 , 1.33212932, 1.54220434, 1.78828965, 2.12824897,
                               #1.25, 1.35, 1.6, 1.80, 2.15
                              ])
    bounds = np.array( 
                       [ [0, 1100], ]   *len(peak_pos_guess) +     #scale
                       [ [1.0,2.6], ]   *len(peak_pos_guess) +     # mean
                       [ [1e-5, 1.0e-1], ] *len(peak_pos_guess) ,     # std
                         ).T
    p0 = np.hstack([np.array([50]*len(peak_pos_guess) ), np.array(peak_pos_guess) , np.array([0.01]*len(peak_pos_guess) )])
    print(p0.shape)
    print(bounds.shape)
    
    # popt= p0
    popt, pcov = curve_fit(gaussians6, 
                           v_prot, 
                           uv_prot, 
                           p0=p0,
                           bounds=bounds ,
                           # maxfev=1e6
                          )
    print(popt)
    peaks = popt[len(peak_pos_guess)+1:len(peak_pos_guess)+len(peak_pos_guess)]
    print(peaks)
    

    # Plot
    fig, ax = plt.subplots(ncols=2, figsize=(12,5))
    ax[0].plot(v_dex, uv_dex, label='Dextran')
    ax[0].plot(v_prot, uv_prot, 'k-', label='Proteins')
    ax[0].plot(v_prot, gaussians6(v_prot, *popt), 'r-', linewidth=1, label='Fit')
    
    for i in range(len(peak_pos_guess)):
        color='r'
        if i == 0:
            color='gray'
        ax[0].fill_between(v_prot, gaussian(v_prot, *popt[i::len(peak_pos_guess)]), color=color, alpha=0.2)

        
    
    ax[0].set(xlim=[0.8,2.5])
    ax[0].legend(loc='center left', bbox_to_anchor=(1, 0.5) )#loc='best')
    ax[0].set(xlabel='Retention volume / mL', ylabel='A280 / a.u.',)# title=chroma.split('/')[-2])

    
    # Calibration curve
    Kav = (peaks - void) / (3-void)
    log10mw = np.log10(list(hmw_std.values()))
    slope, intercept, rvalue, pvalue, stderr = linregress(log10mw, Kav)
    xs = np.linspace(log10mw.min(), log10mw.max(), 100)
    
    for peak in peaks:
        ax[0].axvline(x=peak, color='gray', zorder=-10, alpha=0.5)
        ax[0].text(peak, 1.2*np.max(uv_prot), f"{peak:2.2f} mL", ha="left", va="center", rotation=60, fontsize=10)
    
    ax[0].axvline(x=void, color='k', zorder=-10, alpha=0.8)
    ax[0].text(void, 1.2*np.max(uv_prot), f"Void\n{void:2.2f} mL", ha="right", va="center", rotation=60, fontsize=10)
    
    ax[0].set_ylim( [-1.0, 1.5 * np.max(uv_prot) ] )


    ax[1].scatter(log10mw, Kav, color='white', edgecolor='k')
    ax[1].plot(xs, intercept + slope*xs, 'k-', linewidth=1, label=f'$R^2$ = {rvalue**2:.3f}')
    ax[1].legend(loc='best' ) #loc='best',
    ax[1].set(xlabel='log10(MW)', ylabel='Kav',)# title=chroma.split('/')[-2])
    ax.set_title("S75 ")
    plt.tight_layout()
    plt.savefig(savename + '.png', dpi=300)
    plt.show()
    
    # Save calibration data.
    calibration_data = {
        'log10mw':log10mw.tolist(),
        'Kav':Kav.tolist(),
        'Vc':3,
        'Vo':void,
        'intercept':intercept,
        'slope':slope
    }
    
    print("saving to: " , os.path.abspath ( savename + f'_excl_void_agilent.json' )  )
    with open(savename + f'_excl_void_S200mini_agilent.json' , 'w') as f:
        json.dump(calibration_data, f)

def S200_5_150_calibration_agilent(path_dextran, path_calib, flowrate, comment="" ):
    '''
    chroma: folder path containing the blue-dextran and protein injections
    '''
    
    savename = 'calibrations/' + "HPLC_" + path_dextran.split('/')[-1].replace(".CSV", "").replace(" ", "_") + "_" + comment
    
    hmw_std = {
        'Thyroglobulin':669000,
        'Ferritin':440000,
        'Aldolase':158000,
        'Conalbumin':75000,
        'Ovalbumin':43000,
    }


    dextran = np.loadtxt(path_dextran, delimiter=",")
    dextran[:,0] = dextran[:,0] * flowrate

    v_dex, uv_dex = dextran[:,0] , dextran[:,1]
    void = v_dex[np.argmax(uv_dex)]


    
    print ( f"VOID: {void:1.3f} mL")

    proteins = np.loadtxt(path_calib, delimiter=",")
    proteins[:,0] = proteins[:,0] * flowrate
    v_prot, uv_prot = proteins[:,0], proteins[:,1]

    # plt.plot(dextran[:,0], dextran[:,1])
    # plt.plot(proteins[:,0], proteins[:,1])
    # plt.show()

    # Get peaks by fitting to Gaussian functions
    peak_pos_guess = np.array([#void*0.95, #void peak 
                               1.22892861 , 1.33212932, 1.54220434, 1.7, 1.9,
                               #1.25, 1.35, 1.6, 1.80, 2.15
                              ])
    bounds = np.array( 
                       [ [0, 3200], ]   *len(peak_pos_guess) +     #scale
                       [ [1.0,2.5], ]   *len(peak_pos_guess) +     # mean
                       [ [1e-5, 2.0e-1], ] *len(peak_pos_guess) ,   # std
                         ).T

    p0 = np.hstack([np.array([500]*len(peak_pos_guess) ), np.array(peak_pos_guess) , np.array([0.01]*len(peak_pos_guess) )])
    print(p0.shape)
    print(bounds.shape)
    
    # Plot
    fig, ax = plt.subplots(ncols=2, figsize=(12,5))
    ax[0].plot(v_dex, uv_dex, label='Dextran')
    ax[0].plot(v_prot, uv_prot, 'k-', label='Proteins')
    

    # popt= p0
    popt, pcov = curve_fit(gaussians5, 
                           v_prot, 
                           uv_prot, 
                           p0=p0,
                           bounds=bounds ,
                           # maxfev=1e6
                          )
    print(popt)
    peaks = popt[len(peak_pos_guess):len(peak_pos_guess)+len(peak_pos_guess)]
    print(peaks)
    
    ax[0].plot(v_prot, gaussians5(v_prot, *popt), 'r-', linewidth=1, label='Fit')



    for i in range(len(peak_pos_guess)):
        color='r'
        # if i == 0:
        #     color='gray'
        ax[0].fill_between(v_prot, gaussian(v_prot, *popt[i::len(peak_pos_guess)]), color=color, alpha=0.2)

        
    
    ax[0].set(xlim=[0.8,2.5])
    ax[0].legend(loc='center left', bbox_to_anchor=(1, 0.5) )#loc='best')
    ax[0].set(xlabel='Retention volume / mL', ylabel='A280 / a.u.',)# title=chroma.split('/')[-2])

    
    # Calibration curve
    Kav = (peaks - void) / (3-void)
    log10mw = np.log10(list(hmw_std.values()))
    slope, intercept, rvalue, pvalue, stderr = linregress(log10mw, Kav)
    xs = np.linspace(log10mw.min(), log10mw.max(), 100)
    
    for peak in peaks:
        ax[0].axvline(x=peak, color='gray', zorder=-10, alpha=0.5)
        ax[0].text(peak, 1.2*np.max(uv_prot), f"{peak:2.2f} mL", ha="left", va="center", rotation=60, fontsize=10)
    
    ax[0].axvline(x=void, color='k', zorder=-10, alpha=0.8)
    ax[0].text(void, 1.2*np.max(uv_prot), f"Void\n{void:2.2f} mL", ha="right", va="center", rotation=60, fontsize=10)
    
    ax[0].set_ylim( [-1.0, 1.2 * np.max(uv_prot) ] )


    ax[1].scatter(log10mw, Kav, color='white', edgecolor='k')
    ax[1].plot(xs, intercept + slope*xs, 'k-', linewidth=1, label=f'$R^2$ = {rvalue**2:.3f}')
    ax[1].legend(loc='best' ) #loc='best',
    ax[1].set(xlabel='log10(MW)', ylabel='Kav',)# title=chroma.split('/')[-2])
    plt.tight_layout()
    plt.savefig(savename + '.png', dpi=300)
    print(savename + '.png')
    plt.show()
    
    # Save calibration data.
    calibration_data = {
        'log10mw':log10mw.tolist(),
        'Kav':Kav.tolist(),
        'Vc':3,
        'Vo':void,
        'intercept':intercept,
        'slope':slope
    }
    
    print("saving to: " , os.path.abspath ( savename + f'_S200mini_agilent.json' )  )
    with open(savename + f'_S200mini_agilent.json' , 'w') as f:
        json.dump(calibration_data, f)


#for akta

def S75_5_150_calibration_fitvoid_akta(chroma, no_dextran, no_calib ):
    '''
    chroma: folder path containing the blue-dextran and protein injections
    '''
    
    hmw_std = {
        'Conalbumin':75000,
        'Ovalbumin':43000,
        'Carbonic anhydrase':29000,
        'Ribonuclease A':13700,
        'Aprotinin':6500
    }

    # Dextran injection
    dextran = pc_uni6(glob.glob(f'{chroma}*{no_dextran}.zip')[0])
    dextran.load()
    dextran.xml_parse()
    dextran.clean_up()
    v_dex, uv_dex = np.array(dextran['UV 1_280']['data']).T
    void = v_dex[np.argmax(uv_dex)]

    
    # Proteins injection
    proteins = pc_uni6(glob.glob(f'{chroma}*{no_calib}.zip')[0])
    proteins.load()
    proteins.xml_parse()
    proteins.clean_up()
    v_prot, uv_prot = np.array(proteins['UV 1_280']['data']).T

    # Get peaks by fitting to Gaussian functions
    peak_pos_guess = np.array([1.1, #void peak 
                               1.25, 1.35, 1.64, 1.80, 2.0])
    bounds = np.array( 
                       [ [0,100], ] *len(peak_pos_guess) +
                       [ [1.1,3], ] *len(peak_pos_guess) +
                       [ [1e-3, 1e-1], ] *len(peak_pos_guess) ,
                         ).T
    p0 = np.hstack([np.array([50]*len(peak_pos_guess) ), np.array(peak_pos_guess) , np.array([0.01]*len(peak_pos_guess) )])
    print(p0.shape)
    print(bounds.shape)
    
    popt, pcov = curve_fit(gaussians6, 
                           v_prot, 
                           uv_prot, 
                           p0=p0,
                           bounds=bounds ,
                          )
    #remove void peak
    # peaks = popt[5:10]
    print(popt)
    peaks = popt[len(peak_pos_guess)+1:len(peak_pos_guess)+len(peak_pos_guess)]
    print(peaks)

    # Plot
    fig, ax = plt.subplots(ncols=2, figsize=(12,5))
    ax[0].plot(v_dex, uv_dex, label='Dextran')
    ax[0].plot(v_prot, uv_prot, 'k-', label='Proteins')
    ax[0].plot(v_prot, gaussians6(v_prot, *popt), 'r-', linewidth=1, label='Fit')
    
    for i in range(len(peak_pos_guess)):
        color='r'
        if i == 0:
            color='gray'
        ax[0].fill_between(v_prot, gaussian(v_prot, *popt[i::len(peak_pos_guess)]), color=color, alpha=0.2)

        
    
    ax[0].set(xlim=[0.8,2.5])
    ax[0].legend(loc='best')
    ax[0].set(xlabel='Retention volume / mL', ylabel='A280 / a.u.', title=chroma.split('/')[-2])

    
    # Calibration curve
    Kav = (peaks - void) / (3-void)
    log10mw = np.log10(list(hmw_std.values()))
    slope, intercept, rvalue, pvalue, stderr = linregress(log10mw, Kav)
    xs = np.linspace(log10mw.min(), log10mw.max(), 100)
    
    for peak in peaks:
        ax[0].axvline(x=peak, color='gray', zorder=-10, alpha=0.5)
        ax[0].text(peak, 1.2*np.max(uv_prot), f"{peak:2.2f} mL", ha="left", va="center", rotation=60, fontsize=10)
    
    ax[1].scatter(log10mw, Kav, color='white', edgecolor='k')
    ax[1].plot(xs, intercept + slope*xs, 'k-', linewidth=1, label=f'$R^2$ = {rvalue**2:.3f}')
    ax[1].legend(loc='best')
    ax[1].set(xlabel='log10(MW)', ylabel='Kav', title=chroma.split('/')[-2])
    plt.tight_layout()
    plt.savefig('calibrations/' + chroma.split('/')[-2] + f'_excl_void_{no_calib}.png', dpi=300)
    plt.show()
    
    # Save calibration data.
    calibration_data = {
        'log10mw':log10mw.tolist(),
        'Kav':Kav.tolist(),
        'Vc':3,
        'Vo':void,
        'intercept':intercept,
        'slope':slope
    }
    
    with open('calibrations/' + chroma.split('/')[-2] + f'_excl_void_{no_calib}.json', 'w') as f:
        json.dump(calibration_data, f)

def S200_5_150_calibration_fitvoid(chroma, no_dextran, no_calib ):
    '''
    chroma: folder path containing the blue-dextran and protein injections
    '''
    
    hmw_std = {
        'Thyroglobulin':669000,
        'Ferritin':440000,
        'Aldolase':158000,
        'Conalbumin':75000,
        'Ovalbumin':43000,
    }

    # Dextran injection
    dextran = pc_uni6(glob.glob(f'{chroma}*{no_dextran}.zip')[0])
    dextran.load()
    dextran.xml_parse()
    dextran.clean_up()
    v_dex, uv_dex = np.array(dextran['UV 1_280']['data']).T
    void = v_dex[np.argmax(uv_dex)]

    
    # Proteins injection
    proteins = pc_uni6(glob.glob(f'{chroma}*{no_calib}.zip')[0])
    proteins.load()
    proteins.xml_parse()
    proteins.clean_up()
    v_prot, uv_prot = np.array(proteins['UV 1_280']['data']).T

    # Get peaks by fitting to Gaussian functions
    peak_pos_guess = np.array([1.1, #void peak 
                               1.25, 1.35, 1.64, 1.80, 2.0])
    bounds = np.array( 
                       [ [0,100], ] *len(peak_pos_guess) +
                       [ [1.1,3], ] *len(peak_pos_guess) +
                       [ [1e-3, 1e-1], ] *len(peak_pos_guess) ,
                         ).T
    p0 = np.hstack([np.array([50]*len(peak_pos_guess) ), np.array(peak_pos_guess) , np.array([0.01]*len(peak_pos_guess) )])
    print(p0.shape)
    print(bounds.shape)
    
    popt, pcov = curve_fit(gaussians6, 
                           v_prot, 
                           uv_prot, 
                           p0=p0,
                           bounds=bounds ,
                          )
    #remove void peak
    # peaks = popt[5:10]
    print(popt)
    peaks = popt[len(peak_pos_guess)+1:len(peak_pos_guess)+len(peak_pos_guess)]
    print(peaks)

    # Plot
    fig, ax = plt.subplots(ncols=2, figsize=(12,5))
    ax[0].plot(v_dex, uv_dex, label='Dextran')
    ax[0].plot(v_prot, uv_prot, 'k-', label='Proteins')
    ax[0].plot(v_prot, gaussians6(v_prot, *popt), 'r-', linewidth=1, label='Fit')
    
    for i in range(len(peak_pos_guess)):
        color='r'
        if i == 0:
            color='gray'
        ax[0].fill_between(v_prot, gaussian(v_prot, *popt[i::len(peak_pos_guess)]), color=color, alpha=0.2)

        
    
    ax[0].set(xlim=[0.8,2.5])
    ax[0].legend(loc='best')
    ax[0].set(xlabel='Retention volume / mL', ylabel='A280 / a.u.', title=chroma.split('/')[-2])

    
    # Calibration curve
    Kav = (peaks - void) / (3-void)
    log10mw = np.log10(list(hmw_std.values()))
    slope, intercept, rvalue, pvalue, stderr = linregress(log10mw, Kav)
    xs = np.linspace(log10mw.min(), log10mw.max(), 100)
    
    for peak in peaks:
        ax[0].axvline(x=peak, color='gray', zorder=-10, alpha=0.5)
        ax[0].text(peak, 1.2*np.max(uv_prot), f"{peak:2.2f} mL", ha="left", va="center", rotation=60, fontsize=10)
    
    ax[1].scatter(log10mw, Kav, color='white', edgecolor='k')
    ax[1].plot(xs, intercept + slope*xs, 'k-', linewidth=1, label=f'$R^2$ = {rvalue**2:.3f}')
    ax[1].legend(loc='best')
    ax[1].set(xlabel='log10(MW)', ylabel='Kav', title=chroma.split('/')[-2])
    plt.tight_layout()
    plt.savefig('calibrations/' + chroma.split('/')[-2] + f'_excl_void_{no_calib}.png', dpi=300)
    plt.show()
    
    # Save calibration data.
    calibration_data = {
        'log10mw':log10mw.tolist(),
        'Kav':Kav.tolist(),
        'Vc':3,
        'Vo':void,
        'intercept':intercept,
        'slope':slope
    }
    
    with open('calibrations/' + chroma.split('/')[-2] + f'_excl_void_S200mini_{no_calib}.json', 'w') as f:
        json.dump(calibration_data, f)



