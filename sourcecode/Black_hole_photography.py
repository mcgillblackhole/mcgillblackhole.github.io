import warnings
warnings.filterwarnings('ignore')
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import sys
import matplotlib
matplotlib.use('Qt5Agg')
from PyQt5 import QtCore, QtWidgets
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import scipy as sp
from scipy.integrate import solve_ivp, odeint
import seaborn as sns
from scipy.interpolate import interp1d
from PyQt5.QtCore import QFile, QTextStream
from scipy.special import ellipj
from scipy.special import ellipkinc
from scipy.special import ellipk
from scipy.optimize import fsolve
import time

# ------------------------ Raytracing of light -------------

def S_null(Z, t, p):
    r, rdot, phi = Z
    M, L = p
    phidot = L / r**2
    return [
        rdot,
        L**2 * (r - 3*M) / r**4,
        phidot
    ]


def init_cond(b, x_init):
    r_init = np.sqrt(b**2 + x_init**2)
    phi_init = np.arccos(x_init / r_init)
    rdot_init = np.cos(phi_init)
    phidot_init = -np.sqrt((1 - rdot_init**2) / r_init**2)
    L = r_init**2 * phidot_init
    return [r_init, rdot_init, phi_init], L


def integrate(t, initial, p):
    sol = odeint(S_null, initial, t, args=(p,),)
    r = sol[:,0]
    phi = sol[:,2]
    return r, phi

# ------------------------ Isoradial curves of accretion disk -------------

def B_fun(p, M):
    return (p**3/(p-2*M))**0.5

def Q_fun(P, M):    
    return ((P - 2*M)*(P+6*M))**0.5


def gamma(alpha, theta_0):
    a = np.cos(alpha)**2 + 1/(np.tan(theta_0)**(2))
    return np.arccos( np.cos(alpha)/(a**0.5)  )


def k2(P, M):
    Q = Q_fun(P, M)
    return ((Q-P+6*M)/(2*Q))

def zeta_inf(P, M):
    Q = Q_fun(P, M)
    ratio = (Q-P+2*M)/(Q-P +6*M)
    return np.arcsin( ratio**0.5 )

### 1st order
def Up(P,alpha, M, theta_0): #u = 1/r
    # b = np.sqrt(X**2 + Y**2)
    # alpha = np.arctan(Y/X)
    Q = Q_fun(P, M)
    A1 = (Q-P+2*M)/(4*M*P)
    A2 = (Q-P+6*M)/(4*M*P)
    g = gamma(alpha,theta_0)
    ratio = (Q/P)**0.5
    mod = k2(P, M)
    sn, cn, dn, ph = ellipj( ((g/2)*ratio + ellipkinc(zeta_inf(P, M) , mod)), mod )
    
    return -A1 + A2*sn**2

### 2nd order
def Up2(P,alpha,M,theta_0): #u = 1/r
    # b = np.sqrt(X**2 + Y**2)
    # alpha = np.arctan(Y/X)
    Q = Q_fun(P, M)
    A1 = (Q-P+2*M)/(4*M*P)
    A2 = (Q-P+6*M)/(4*M*P)
    g = gamma(alpha,theta_0)
    ratio = (Q/P)**0.5
    mod = k2(P,M)
    sn, cn, dn, ph = ellipj( (0.5*(g-2*np.pi)*ratio + 2*ellipk(mod) - ellipkinc(zeta_inf(P, M) , mod)), mod )
    
    return -A1 + A2*sn**2
    
def plot_solution(self, fig, x_1, y_1, xs_acc_disk, ys_acc_disk, xs_acc_disk_2nd, ys_acc_disk_2nd, angle_disk, accret_disk_min_r, accret_disk_max_r, M, t):
    
    ax1 = fig.add_subplot(1,2,1)
    #ax1.figure.clear()
    ax1.grid(False)
    x_disk_plus = np.arange(accret_disk_min_r, accret_disk_max_r, 0.05)
    x_disk_minus = np.arange(-accret_disk_max_r, -accret_disk_min_r, 0.05)
    for i in range(len(x_1)):
        attempts = 0
        for j in range(len(x_disk_plus)):
            attempts += 1
            d = np.sqrt((x_1[i] - x_disk_plus[j])**2 + (y_1[i]-x_disk_plus[j]*np.tan(angle_disk))**2)
            r = np.sqrt(x_1[i]**2 + y_1[i]**2)
            max_idx = np.where(r < 2*M)[0]
            if np.min(d) < 0.1:
                idx = j
                plt.scatter(x_disk_plus[idx], x_disk_plus[idx]*np.tan(angle_disk), color = 'black', s = 20)
                line_hit, = plt.plot(x_1[i], y_1[i], color = 'blue')
                break
        
        attempts_2 = 0
        for j in range(len(x_disk_minus)):
            attempts_2 += 1
            d = np.sqrt((x_1[i] - x_disk_minus[j])**2 + (y_1[i]-x_disk_minus[j]*np.tan(angle_disk))**2)
            if np.min(d) < 0.1:
                idx = j
                plt.scatter(x_disk_minus[idx], x_disk_minus[idx]*np.tan(angle_disk), color = 'black', s = 20)
                plt.plot(x_1[i], y_1[i], color = 'blue')
                break

        if attempts == len(x_disk_plus) and attempts_2 == len(x_disk_minus):
            line_fail, = plt.plot(x_1[i], y_1[i], color = 'lightgrey')
            # time.sleep(1)
            # plt.close()
    circle = plt.Circle((0., 0.), 2*M, color='black', fill=True, zorder=10)
    # plt.legend()
    #axes.add_artist(circle)
    plt.gca().add_patch(circle)
    # plt.plot(x_disk,[0 for i in range(len(x_disk))], color = 'black', linestyle = 'dashed')
    plt.plot(x_disk_plus, x_disk_plus*np.tan(angle_disk), color = 'red')
    plt.plot(x_disk_minus, x_disk_minus*np.tan(angle_disk), color = 'red')
    plt.xlabel(r"Distance $X$", fontsize = 10, color = 'white')
    plt.ylabel(r"Distance $Z$", fontsize = 10, color = 'white')
    plt.yticks(color = 'white')
    plt.xticks(color = 'white')
    plt.axis('equal')
    # legend
    try:
        ax1.legend(handles = [line_hit, line_fail], labels = ['Accretion disk hit', 'Off to infinity or captured'], loc='upper right')
    except UnboundLocalError:
        print('All lightrays either hit the accretion disk or miss it.')

    ax1.set_xlim(-20,20)
    ax1.set_xbound(-20,20)
    ax1.set_ylim(-20,20)
    ax1.set_ybound(-20,20)
    # set title
    ax1.set_title(r"Lightlike geodesics for" "\n" r"$ds^2 = -(1-2M/r) \ dt^2 + (1 - 2M/r)^{-1}dr^2 + r^2 d\Omega^2$",
                fontsize=12, color = 'white')
        
    # plot accretion disk view
    ax2 = fig.add_subplot(1,2,2)
    #ax1.figure.clear()
    ax2.grid(False)

    for i in range(len(xs_acc_disk)):
        x = xs_acc_disk[i]
        y = ys_acc_disk[i]
        x2 = xs_acc_disk_2nd[i]
        y2 = ys_acc_disk_2nd[i]
        if np.pi/2 - angle_disk < 78*np.pi/180:
            if i > 3:
                x_interp = np.concatenate((y[-20:-10], y[10:20]))
                y_interp = np.concatenate((-x[-20:-10], -x[10:20]))
                f = interp1d(x_interp, y_interp, kind = 'quadratic')

                x_new = np.linspace(y[-10], y[10], 100)
                plt.plot(x_new, f(x_new), color = sns.color_palette('hot_r', len(xs_acc_disk))[i],
                            zorder = 2, lw = 4)
    
        plt.plot(y2, x2, color = sns.color_palette('hot_r', len(xs_acc_disk))[i],
                zorder = 1, lw = 4)
        plt.plot(y,-x, color = sns.color_palette('hot_r', len(xs_acc_disk))[i],
                zorder = 2, lw = 4)

    circle = plt.Circle((0., 0.), 2*M, color='black', fill=True, zorder=1)
    plt.gca().add_patch(circle)
    circle = plt.Circle((0., 0.), 3*M, color='black', fill=True, zorder=1, alpha = 0.3)
    plt.gca().add_patch(circle)
    circle = plt.Circle((0., 0.), np.sqrt(3)*3*M, facecolor = 'none', edgecolor='white', fill=True,
                        zorder=1, lw = 4)
    plt.gca().add_patch(circle)

    ax2.set_xlim(-32,32)
    ax2.set_ylim(-32,32)
    ax2.axis('off')
    fig.patch.set_color('dimgrey')
    for ax in fig.axes:
        ax.patch.set_color('dimgrey')
    ax2.set_title(f'Black hole accretion disk view at {int((np.pi/2 - angle_disk)*180/np.pi)}º', fontsize = 14,
            color = 'white')

    print('Done! Check your plots :)')

# --------------------------------------------------------------------------------- #

class Widget(QtWidgets.QWidget):
    def __init__(self, *args, **kwargs):
        QtWidgets.QWidget.__init__(self, *args, **kwargs)
        global fig
        fig = plt.figure(figsize=(3, 3), dpi=100)    
        fig.set_facecolor('#021226')
        self.figure = fig
        self.canvas = FigureCanvas(self.figure)
        self.canvas.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)

        #button for starting simulation 
        button = QtWidgets.QPushButton("Run Simulation!")
        button.clicked.connect(self.main)
        #button for stopping simulation 
        button2 = QtWidgets.QPushButton("Stop Simulation")
        button2.clicked.connect(self.stop)
        #button for starting simulation 
        button = QtWidgets.QPushButton("Random Simulation!")
        button.clicked.connect(self.main)
        #button for stopping simulation 
        button2 = QtWidgets.QPushButton("Stop Simulation")
        button2.clicked.connect(self.stop)
        #button for choosing initial conditions
        combobox = QtWidgets.QComboBox(self)
        # Add preconfigured cases to the box
        combobox.addItem("Loop")
        choose = QtWidgets.QLabel("Choose this predefined configuration to see what happens to a lightray when the impact parameter is critic! (Don't forget to stop the simulation before trying another config).")
        combobox.activated[str].connect(self.case) 
        #button.resize(150, 50)
        lay = QtWidgets.QVBoxLayout(self)
        vlay1 = QtWidgets.QHBoxLayout(self)
        vlay2 = QtWidgets.QHBoxLayout(self)
        vlay3 = QtWidgets.QHBoxLayout(self)
        vlay4 = QtWidgets.QHBoxLayout(self)
        vlay5 = QtWidgets.QHBoxLayout(self)
        vlay6 = QtWidgets.QHBoxLayout(self)
        #lay = QtWidgets.QGridLayout(self)
        lay.addWidget(self.canvas)
        #vlay1.addWidget(button)

        #vlay2.addWidget(button3)
        #vlay2.addWidget(button4)
        vlay1.addWidget(choose)
        vlay2.addWidget(combobox)
        vlay2.addWidget(button2)
        groupBox = QtWidgets.QGroupBox()
        groupBox.setTitle('McGill Physics Hackathon 2022: Black Hole Photography')
        groupBox.setStyleSheet(self.getStyleSheet("./styles.css"))
        groupBox.setAlignment(QtCore.Qt.AlignCenter)
        lay.addLayout(vlay1)
        lay.addLayout(vlay2)
        lay.addLayout(vlay3)
        lay.addLayout(vlay4)
        lay.addLayout(vlay5)
        lay.addLayout(vlay6)
        groupBox.setLayout(lay)

        choose2 = QtWidgets.QLabel("You can also choose the mass, impact parameters and the angle between the observer and the accretion disk, and run a custom built simulation.")
        
        # Sliders for the masses
        a_slider = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        a_slider.setRange(0,50)
        a_slider.setValue(1)
        self.a_slider = a_slider

        button = QtWidgets.QPushButton("Ready, Set, Go!")
        button.clicked.connect(self.main)

        vlay3.addWidget(choose2)
        vlay5.addWidget(button)
        vlay5.addWidget(button2)
        b_slider = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        b_slider.setRange(0,50)
        b_slider.setValue(1)
        self.b_slider = b_slider

        c_slider = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        c_slider.setRange(0,50)
        c_slider.setValue(1)
        self.c_slider = c_slider

        a_value= QtWidgets.QLabel()
        b_value = QtWidgets.QLabel() 
        c_value = QtWidgets.QLabel() 

        main_layout = QtWidgets.QGridLayout(self)
        main_layout.addWidget(groupBox,0,0)
        a_slider.valueChanged.connect(a_value.setNum)
        b_slider.valueChanged.connect(b_value.setNum)
        c_slider.valueChanged.connect(c_value.setNum)   

        # Values for masses
        self.spinBoxm8 = QtWidgets.QDoubleSpinBox()
        self.spinBoxm8.setRange(1, 10)
        self.spinBoxm8.setValue(1)
        self.spinBoxm8.setPrefix("BH mass =  ")
        self.spinBoxm8.valueChanged.connect(c_value.setNum)
        vlay4.addWidget(self.spinBoxm8)

        # Values b #vamos ter que mudar esses abc_value lÃ¡ em cima
        self.spinBoxm1 = QtWidgets.QDoubleSpinBox()
        self.spinBoxm1.setRange(-20, 20)
        self.spinBoxm1.setValue(3*np.sqrt(3)*self.spinBoxm8.value())
        self.spinBoxm1.setPrefix("Minimum impact parameter =  ")
        self.spinBoxm1.valueChanged.connect(a_value.setNum)
        vlay4.addWidget(self.spinBoxm1)

        self.spinBoxm2 = QtWidgets.QDoubleSpinBox()
        self.spinBoxm2.setRange(-20, 20)
        self.spinBoxm2.setValue(3*np.sqrt(3)*self.spinBoxm8.value())
        self.spinBoxm2.setPrefix("Maximum impact parameter =  ")
        self.spinBoxm2.valueChanged.connect(b_value.setNum)
        vlay4.addWidget(self.spinBoxm2)

        self.spinBoxm3 = QtWidgets.QDoubleSpinBox()
        self.spinBoxm3.setRange(0, 50)
        self.spinBoxm3.setValue(1)
        self.spinBoxm3.setPrefix("Number of lightrays =  ")
        self.spinBoxm3.valueChanged.connect(c_value.setNum)
        vlay4.addWidget(self.spinBoxm3)

        self.spinBoxm4 = QtWidgets.QDoubleSpinBox()
        self.spinBoxm4.setRange(-40, 40)
        self.spinBoxm4.setValue(-40)
        self.spinBoxm4.setPrefix("Initial position =  ")
        self.spinBoxm4.valueChanged.connect(c_value.setNum)
        vlay4.addWidget(self.spinBoxm4)

        self.spinBoxm7 = QtWidgets.QDoubleSpinBox()
        self.spinBoxm7.setRange(0, 90)
        self.spinBoxm7.setValue(10)
        self.spinBoxm7.setPrefix("Theta =  ")
        self.spinBoxm7.valueChanged.connect(c_value.setNum)
        vlay4.addWidget(self.spinBoxm7)

        #lay.addWidget(button2)
        #self.plot()


    def case(self, text):
        if text == "Loop":
            self.Loop()

    global anispeed 
    anispeed = 0.0002

    
    def getStyleSheet(self, path):
        f = QFile(path)
        f.open(QFile.ReadOnly | QFile.Text)
        stylesheet = QTextStream(f).readAll()
        f.close()
        return stylesheet

    def stop(self):
        ani.event_source.stop()
        ani.frame_seq = ani.new_frame_seq()  

    def Loop(self):
        #ani.event_source.stop()
        self.figure.clear()
        M = 1
        angle_disk = 0
        theta_0 = np.pi/2
        accret_disk_min_r = 6*M
        accret_disk_max_r = 15*M
        r_list = np.arange(4,30,0.5)*M
        b_c = 3*np.sqrt(3)*M
        xs_acc_disk = []
        ys_acc_disk = []

        for idx, r in enumerate(r_list):
            alpha_list = np.arange(0, 4*np.pi,0.01)
            b_list = []
            alpha_res = []
            #plt.figure(figsize = (15,15), dpi = 300)
            for alpha in alpha_list:    
                def cu(P,M,theta_0):
                    return 1 - r*Up(P,alpha,M,theta_0)
                
                #P_list = np.arange(3*M, 50*M,0.1)
                
                #plt.plot(P_list, cu(P_list))
                
                guess = b_c + 0.1
                raiz = fsolve(cu,[guess], args = (M,theta_0))[0]
                if raiz != guess:
                
                    #plt.scatter(raiz, cu(raiz))
                    
                    b_raiz = B_fun(raiz, M)
                    b_list.append(b_raiz)
                    alpha_res.append(alpha)
            
            x = b_list*np.cos(alpha_res)
            y = b_list*np.sin(alpha_res)
            xs_acc_disk.append(x)
            ys_acc_disk.append(y)
        # time interval and the step size
        x_init = -50
        max_t = 100.
        # impact parameters to plot
        b_min = 3*np.sqrt(3)*M
        b_max = 3*np.sqrt(3)*M
        b_num = 1
        b_list = np.linspace(b_min, b_max, int(b_num))
        # time points to evaluate
        t = np.arange(0, max_t, 0.01)
        # t = np.arange(0, 30, anispeed)
        
        # vectors for the solutions
        x_1 = np.zeros((len(t)))
        y_1 = np.zeros((len(t)))

        # initial conditions
        for j, b in enumerate(b_list):
            initial, L = init_cond(b, x_init)
    
            r, phi = integrate(t, initial, [M, L])

            if b < 0:
                x_1 = r * np.cos(phi)
                y_1 = - r * np.sin(phi)
            else:
                x_1 = r * np.cos(phi)
                y_1 = r * np.sin(phi)
        
        plot_solution(self, fig, [x_1], [y_1], xs_acc_disk, ys_acc_disk, angle_disk, accret_disk_min_r, accret_disk_max_r, M, t)
        
    def main(self):
        self.figure.clear()
        M = self.spinBoxm8.value()
        angle_disk = self.spinBoxm7.value()*np.pi/180
        theta_0 = np.pi/2 - self.spinBoxm7.value()*np.pi/180
        accret_disk_min_r = 6*M
        accret_disk_max_r = 15*M
        r_list = np.arange(4,30,0.5)*M
        b_c = 3*np.sqrt(3)*M
        xs_acc_disk = []
        ys_acc_disk = []
        xs_acc_disk_2nd = []
        ys_acc_disk_2nd = []

        for idx, r in enumerate(r_list):
            if theta_0 < 78*np.pi/180:
                alpha_list = np.arange(0, 2*np.pi,0.01)
            else:
                alpha_list = np.arange(0, 4*np.pi,0.01)
            b_list = []
            b_list2 = []
            alpha_res = []
            alpha_res2 = []
            #plt.figure(figsize = (15,15), dpi = 300)
            for alpha in alpha_list:    
                # 1st order figure
                def cu(P,M,theta_0):
                    return 1 - r*Up(P,alpha,M,theta_0)

                # 2nd order figure
                def cu2(P,M,theta_0):
                    return -(1 - r*Up2(P,alpha,M,theta_0))


                guess = b_c + 0.1
                # raiz = fsolve(cu,[guess])[0]
                raiz = fsolve(cu,[guess], args = (M,theta_0))[0]
                raiz2 = fsolve(cu2,[guess], args = (M,theta_0))[0]

                # 1st order
                if raiz != guess:

                    b_raiz = B_fun(raiz, M)
                    b_list.append(b_raiz)
                    alpha_res.append(alpha)

                # 2nd order
                if raiz2 != guess:

                    b_raiz2 = B_fun(raiz2, M)
                    b_list2.append(b_raiz2)
                    alpha_res2.append(alpha)

            # 1st order:
            x = b_list*np.cos(alpha_res)
            y = b_list*np.sin(alpha_res)
            xs_acc_disk.append(x)
            ys_acc_disk.append(y)


            #2nd order:
            x2 = b_list2*np.cos(alpha_res2)
            y2 = b_list2*np.sin(alpha_res2)

            ## Limits
            xlims = [-32,32]
            ylims = [-32,32]

            xplot = []
            yplot = []

            for xp in x2:
                if xlims[0] < xp < xlims[1]:
                    xplot.append(xp)

            for yp in y2:
                if ylims[0] < yp < ylims[1]:
                    yplot.append(yp)

        ## Remove points
            if abs(np.shape(xplot)[0] - np.shape(yplot)[0]) != 0:

                xshape = np.shape(xplot)[0] 
                yshape = np.shape(yplot)[0]

                if xshape < yshape:
                    dif = abs(np.shape(xplot)[0] - np.shape(yplot)[0]) 

                    yplot2 = yplot[dif:]
                    xplot2 = xplot

                if yshape < xshape:
                    dif = abs(np.shape(xplot)[0] - np.shape(yplot)[0]) 

                    xplot2 = xplot[dif:]
                    yplot2 = yplot 

            xs_acc_disk_2nd.append(xplot2)
            ys_acc_disk_2nd.append(yplot2)
            


        # time interval and the step size
        x_init = -50
        max_t = 100.
        # impact parameters to plot
        b_min = self.spinBoxm1.value()
        b_max = self.spinBoxm2.value()
        b_num = self.spinBoxm3.value()
        b_list = np.linspace(b_min, b_max, int(b_num))
        # time points to evaluate
        t = np.arange(0, max_t, 0.01)
        # t = np.arange(0, 30, anispeed)
        
        # vectors for the solutions
        x_1 = np.zeros((len(t)))
        y_1 = np.zeros((len(t)))

        # initial conditions
        xs = []
        ys = []
        for j, b in enumerate(b_list):
            initial, L = init_cond(b, x_init)
    
            r, phi = integrate(t, initial, [M, L])

            if b < 0:
                x_1 = r * np.cos(phi)
                y_1 = - r * np.sin(phi)
            else:
                x_1 = r * np.cos(phi)
                y_1 = r * np.sin(phi)
            xs.append(x_1)
            ys.append(y_1)
        
        plot_solution(self, fig, xs, ys, xs_acc_disk, ys_acc_disk, xs_acc_disk_2nd, ys_acc_disk_2nd, angle_disk, accret_disk_min_r, accret_disk_max_r, M, t)
           
        
if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    w = Widget()
    w.show()
    #main()
    sys.exit(app.exec_())
