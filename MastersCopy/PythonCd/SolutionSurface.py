'''
Set of functions to facricate a 2D space containing combinations of real and 
imaginary parts of the wavenumber k_. (analogously E and eta) Its a brute force
solution instead of using the bobyqa in MBBMs code. 
'''
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import BeamWaves as fW

def GetKcomplex(eta, E, I, freq, mu):
    omega = freq*2*np.pi
    D= E*I
    D_ = D*(1+1j*eta)
    k_= pow(mu * omega * omega / D_,0.25) #me faltaba aqui el cuadrado del omega
    return k_
    
def explore(cb, eta, freq, Xpos, meas):

    k_ = GetKcomplex(eta,)
    k_ = 2*np.pi*freq /cb * (1-1j*eta/4)
    #define matrix A from k complex
    kivect = np.zeros((1,4),dtype = complex)
    kivect[0,0] = -1j*k_
    kivect[0,1] = -k_
    kivect[0,2] = 1j*k_
    kivect[0,3] = k_
    kivect = np.outer(Xpos,kivect)

    A = np.ones((Xpos.size,4),dtype=complex)
    A = A * np.e
    A = A**kivect

    s= meas.size
    b=np.zeros((s,1), dtype=complex)
    b[:,0] = meas[:]
    coeffs, residuals, rank, sing = np.linalg.lstsq(A,b, rcond = None)
    return coeffs, residuals, rank, sing

def explore2(cb, eta, freq, Xpos, meas, mu,I):

    Ecalc = np.power(cb,4)*mu/I/(np.power(freq*np.pi*2,2))
    k_ = GetKcomplex(eta, Ecalc, I,freq, mu)
    #k_ = 2*np.pi*freq /cb * (1-1j*eta/4)
    #define matrix A from k complex
    kivect = np.zeros((1,4),dtype = complex)
    kivect[0,0] = -1j*k_
    kivect[0,1] = -k_
    kivect[0,2] = 1j*k_
    kivect[0,3] = k_
    kivect = np.outer(Xpos,kivect)

    A = np.ones((Xpos.size,4),dtype=complex)
    A = A * np.e
    A = A**kivect
    
    #print("conditioning number of A: "+str(np.linalg.cond(A)))

    s= meas.size
    b=np.zeros((s,1), dtype=complex)
    b[:,0] = meas[:]
    #coeffs, residuals, rank, sing = np.linalg.lstsq(A,b, rcond = None)
    #Normalisierungs Versuch
    column_norms = np.linalg.norm(A,axis=0)#this normalization is new because from olf measurements it wasnt working without it, maybe because of big x?
    Anorm  = A / np.maximum(column_norms, 1e-10)
    #print("conditioning number of Anorm: "+str(np.linalg.cond(Anorm)))
    coeffsnorm, residuals2, rank2, sing2 = np.linalg.lstsq(Anorm,b, rcond = None)
    coeffs2 = coeffsnorm / column_norms
    #residuals2[0] = (np.linalg.norm(b-A@coeffs2)**2) #benutz logarithmische Residuen wegen sehr große Zahlen
    
    return coeffs2, residuals2, rank2, sing2

def openFEMresults(filename):
    # Step 1: Detect the line with column titles
    with open(filename, 'r') as file:
        for i, line in enumerate(file):
            if not line.startswith('#'):
                column_title_line = i
                break

    # Step 2: Load the data using Pandas
    data = pd.read_csv(filename, sep='\s+', skiprows=column_title_line)

    # Extract frequency and responses
    frequencies = data.iloc[:, 0].values  # First column (frequencies)
    responses_real = data.iloc[:, 1::2].values  # Real parts (2nd and 4th columns)
    responses_imag = data.iloc[:, 2::2].values  # Imaginary parts (3rd and 5th columns)

    # Combine real and imaginary parts into a complex matrix
    responses = responses_real + 1j * responses_imag
    return frequencies, responses


def exploreSpace(meas, freq, h,b,rho, rangeCb, rangeEta, resolCb, resoleEta,Xpos, plot = ""):

    I = b*h*h*h/12
    mu = h*b*rho

    #Define X position vector
    # start = 0.170
    # step = 0.020
    # x_size = 13
    # stop = start + (x_size * step)
    # Xpos = np.arange(start,stop,step)


    # Define the range for cb values
    cb_n = resolCb
    cb_values = np.linspace(rangeCb[0], rangeCb[1], cb_n)

    # Define the range of eta values
    eta_n = resoleEta
    eta_values = np.linspace(rangeEta[0], rangeEta[1], eta_n)  # Avoid 0 to prevent division errors

    #Mesh grid
    X = cb_values
    Y = eta_values
    X,Y = np.meshgrid(X,Y, indexing = 'ij')

    # Solution space matrix
    Sol = np.zeros((cb_n,eta_n),dtype=float)


    for i in range(cb_values.size):
        for j in range(eta_values.size):

            cb = X[i,j]
            eta = Y[i,j]
            coefs, residuals, rank, sing = explore2(cb, eta, freq, Xpos, meas, mu, I)

            if residuals.size > 0:
                Sol[i,j] = residuals[0]
            else:
                Sol[i,j] = None
    Sol = np.ma.masked_where(np.isnan(Sol), Sol)

    if plot != "":

        #Plot the Error in the solution space


        fig,ax = plt.subplots(subplot_kw={"projection":"3d"})
        surf = ax.plot_surface(X,Y,Sol, vmin = Sol.min(), cmap = cm.viridis)
        # Set labels for each axis
        ax.set_xlabel('cb')
        ax.set_ylabel('eta')
        ax.set_zlabel('Residual')

        ax.set_xticks(np.linspace(X.min(), X.max(), 5))  # Customize the number of ticks
        ax.set_yticks(np.linspace(Y.min(), Y.max(), 5))  # Customize the number of ticks
        ax.set_zticks(np.linspace(Sol.min(), Sol.max(), 5))  # Customize the number of ticks

        fig.colorbar(surf, ax=ax, shrink=0.5, aspect=5, pad=0.1)
        plt.title(plot)
        plt.show()

    min_index = np.unravel_index(np.argmin(Sol),Sol.shape)
    min_cb = cb_values[min_index[0]]
    min_eta = eta_values[min_index[1]]

    coefs, residuals, rank, sing = explore2(min_cb, min_eta, freq, Xpos, meas, mu, I)


    #define beam parameters

    E_ideal = 210000000000
    eta_soll = 0.3
    E_imag = 1j*eta_soll*E_ideal
    E_mag = np.abs(E_ideal + E_imag)
    omega = 2 * np.pi * freq
    rho = 7850
    Ecalc = np.power(min_cb,4)*mu/I/(np.power(omega,2))
    #Ecalc = np.power(min_cb,4)*mu/I/(np.power(freq*np.pi*2,2))
    min_cb_ideal = np.power((E_ideal*I*omega*omega/mu),0.25)

    print("Eta = ", min_eta)
    print("cb = ", min_cb)
    print("cb_ideal = ", min_cb_ideal)
    print("Ecalc = ", Ecalc)

    return min_cb, min_eta, coefs, residuals
if __name__ == "__main__":
    # exploreSpace("res0p800_3mm.txt", 350, 0.003, (115,125), (0.65,0.7), 500,500)
    # print(" ")
    # exploreSpace("res0p800_6mm.txt", 350, 0.006, (165,175), (0.65,0.7), 500,500)
    # print(" ")
    # exploreSpace("res0p290_schm.txt", 350, 0.001, (63,67), (0.27,0.29), 500,500)

    # exploreSpace("res0p100_schm.txt", 350, 0.003, (60,66), (0.05,0.15), 500,500)
    # print(" ")
    # exploreSpace("res0p300_lang.txt", 350, 0.006, (60,66), (0.28,0.31), 500,500)
    # print(" ")
    # exploreSpace("res0p300_lang_1El.txt", 350, 0.001, (60,66), (0.28,0.31), 500,500)
    frequencies, responses = openFEMresults("input/res0p800_1mm_nu0.txt")
    num= 350
    freq = frequencies[num]
    meas = responses[num,:]

    #Define X position vector
    start = 0.070
    step = 0.020
    x_size = 21
    stop = start + (x_size * step)
    Xpos = np.arange(start,stop,step)

    # min_cb, min_eta,coefs = exploreSpace(meas, freq, 0.001, (50,90), (0.6,0.81), 500,500, "0p800")

    # fabricate measurements
    VMs = fW.fabricateVibration(eta=0.0, cb=50, freq=freq, Vp=1, Vpj=0.1, Vm=0.1, Vmj=0.02, plot="Identified from Meas")
    # def fabricateVibration(eta=0.001, cb=65, freq=500, Vp=1, Vm=0, Vpj=0, Vmj=0, plot=""):
    VMs = VMs.flatten()
    min_cb, min_eta, coefs2 = exploreSpace(VMs, freq, 0.001, (10,100), (0.1,1), 500,500, "Re-Identified Error Space")




    # frequencies, responses = openFEMresults("res0p800_1mm_nu0p1.txt")
    # num= 350
    # freq = frequencies[num]
    # meas = responses[num,:]
    # min_cb2, min_eta2 = exploreSpace(meas, freq, 0.001, (50,90), (0.6,0.81), 500,500, "0p800")

    # n = np.arange(300,350,1)
    # etas = np.zeros(51)
    # cbs = np.zeros(51)
    # i=0;
    # for ii in n:
    #     cbi, etai = exploreSpace("res0p290_schm.txt", ii, 0.003, (30,300),(0.15,0.40), 250,250, plot = False)
    #     etas[i] = etai
    #     cbs[i] = cbi
    #     i=i+1

