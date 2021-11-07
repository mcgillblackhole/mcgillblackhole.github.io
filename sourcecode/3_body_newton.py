import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D
import sys
import matplotlib
matplotlib.use('Qt5Agg')
from PyQt5 import QtCore, QtWidgets
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import scipy as sp
import scipy.integrate
from scipy.integrate import solve_ivp
from PyQt5.QtCore import QFile, QTextStream

G=1#6.67408e-11 #N-m2/kg2

eps12 = 1e-8
eps13 = 1e-8
eps23 = 1e-8

# ------------------------ Post Newtonian Solution -------------

def dist(x1, x2):
    
    d = np.sqrt(np.sum((x1- x2)**2)+1e-4) # numerical cutoff of 1e-4 to avoid masses going to infinity
    
    return d

def force(m2, x1, x2):
    
    f0 = m2*(x2 - x1)/dist(x1, x2)**3
    
    return f0

def force2(m1, m2, m3, x1, x2, x3):
    
    dis = m2/dist(x1, x2)
    f2 = force(m1, x2,x1) + force(m3, x2,x3)
    
    return (7/2)*dis*f2

def force3(m2, x1, x2, v1, v2):
    
    f = m2/(dist(x1,x2)**2)*(np.dot((x1 - x2)/dist(x1,x2),(4*v1 - 3*v2)) * (v1-v2))
    
    return f

def force4(m1, m2, m3, x1, x2, x3, v1, v2):
    
    f = force(m2, x1, x2)*(np.dot(v1, v1) + 2*np.dot(v2,v2) - 4*np.dot(v1,v2) - 
                       3/2 * np.dot((x1-x2)/dist(x1,x2), v2)**2 + 
                       1/2 * np.dot(x2 - x1, force(m1, x2, x1) + force(m3, x2, x3)))
    return f

def force5(m2, x1, x2, x3):
    
    r12 = 1/dist(x1, x2)
    r23 = 1/dist(x2, x3)
    r13 = 1/dist(x1, x3)
    f1 = force(m2, x1,x2)*(r12 + r13) 
    f2 = force(m2, x1,x2)*(r12 + r23)
    
    return (-4*f1 -f2)

def dydx(t, u, m1, m2, m3):
    
    """
    x1, y1, z1, vx1, vy1, vz1 = u[:6]
    x2, y2, z2, vx2, vy2, vz2 = u[6:12]
    x3, y3, z3, vx3, vy3, vz3 = u[12:]
    """
    
    x1 = u[:3]
    x2 = u[6:9]
    x3 = u[12:15]
    
    v1 = u[3:6]
    v2 = u[9:12]
    v3 = u[15:]
    
    f01 = force(m2, x1, x2) + force(m3, x1, x3)
    f11 = force2(m1, m2, m3, x1, x2, x3) + force2(m1, m3, m2, x1 ,x3, x2)
    f21 = force3(m2, x1, x2, v1, v2) + force3(m3, x1, x3, v1, v3)
    f31 = force4(m1, m2, m3, x1, x2, x3, v1, v2) + force4(m1, m3, m2, x1, x3, x2, v1, v3)
    f41 = force5(m2, x1,x2,x3) + force5(m3, x1,x3,x2)

    
    f02 = force(m1, x2,x1) + force(m3, x2, x3)
    f12 = force2(m2, m1, m3, x2, x1, x3) + force2(m2, m3, m1, x2, x3, x1)
    f22 = force3(m1, x2, x1, v2, v1) + force3(m3, x2, x3, v2, v3)
    f32 = force4(m2, m1, m3, x2, x1, x3, v2, v1) + force4(m2, m3, m1, x2, x3, x1, v2, v3)
    f42 = force5(m1, x2, x1, x3) + force5(m3, x2, x3, x1)

    
    f03 = force(m1, x3, x1) + force(m2, x3, x2)
    f13 = force2(m3, m1, m2, x3, x1, x2) + force2(m3, m2, m1, x3, x2, x1)
    f23 = force3(m1, x3, x1, v3, v1) + force3(m2, x3, x2, v3, v2)
    f33 = force4(m3, m2, m1, x3, x2, x1, v3, v2) + force4(m3, m1, m2, x3, x1, x2, v3, v1)
    f43 = force5(m1, x3,x1,x2) + force5(m2, x3,x2,x1)
    
    tol2 = 1e-3
    tol3 = 1e-2
    tol4 = 1e-2
    tol4 = 1e-3
    
    return np.array([*v1,*(f01 + tol2*f11 + tol3*f21 + tol4*f31 + tol4*f41 ),  ## 1,2 + 1,3
                     *v2,*(f02 + tol2*f12 + tol3*f22 + tol4*f32 + tol4*f42 ),  ## 2,1 + 2,3
                     *v3,*(f03 + tol2*f13 + tol3*f23 + tol4*f33 + tol4*f43 )])

G = M = 1

# --------------------- Newtonian Solution -----------------------------

def ThreeBodyEquations(t,w,G,m1,m2,m3):
    r1=w[:3]
    r2=w[3:6]
    r3=w[6:9]
    v1=w[9:12]
    v2=w[12:15]
    v3=w[15:18]
    
    r12=sp.linalg.norm(r2-r1)+eps12
    r13=sp.linalg.norm(r3-r1)+eps13
    r23=sp.linalg.norm(r3-r2)+eps23
    
    dv1bydt=m2*(r2-r1)/r12**3+m3*(r3-r1)/r13**3
    dv2bydt=m1*(r1-r2)/r12**3+m3*(r3-r2)/r23**3
    dv3bydt=m1*(r1-r3)/r13**3+m2*(r2-r3)/r23**3
    dr1bydt=v1
    dr2bydt=v2
    dr3bydt=v3
    r12_derivs=np.concatenate((dr1bydt,dr2bydt))
    r_derivs=np.concatenate((r12_derivs,dr3bydt))
    v12_derivs=np.concatenate((dv1bydt,dv2bydt))
    v_derivs=np.concatenate((v12_derivs,dv3bydt))
    derivs=np.concatenate((r_derivs,v_derivs))
    return derivs
    
def plot_solution(self, fig, x_1, y_1, z_1, t, x_2, y_2, z_2, x_3, y_3, z_3, 
         x_4, y_4, z_4, x_5, y_5, z_5, x_6, y_6, z_6):
    
    ax1 = fig.add_subplot(1, 2, 1, projection='3d')
    #ax1.figure.clear()
    ax1.grid(False)
    line, = ax1.plot([], [], [], color='#A325D9', linewidth=0.5, label='$m_{1}$')
    line2, = ax1.plot([], [], [], color='#04AD7B', linewidth=0.5, label='$m_{2}$')
    line3, = ax1.plot([], [], [], color='#D93425', linewidth=0.5, label='$m_{3}$')
    
    point, = ax1.plot([], [], [], marker='o', color='#A325D9', markersize=4)
    point2, = ax1.plot([], [], [], marker='o', color='#04AD7B', markersize=4)
    point3, = ax1.plot([], [], [], marker='o', color='#D93425', markersize=4)
    
    ax1.set_facecolor('#021226')
    ax1.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.2))
    ax1.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.2))
    ax1.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.2))
    # Take out the axis
    #ax1.set_axis_off()
    
    # set point-of-view: specified by (altitude degrees, azimuth degree)
    ax1.view_init(elev=23, azim=-131)
    
    # set limits 
    ax1.set_xlim(-2, 2)
    ax1.set_ylim(-2, 2)
    ax1.set_zlim(-2, 2)
    
    # # set axis labels
    # ax1.set_xlabel('X', fontsize=10)
    # ax1.set_ylabel('Y', fontsize=10)
    # ax1.set_zlabel('Z', fontsize=10)
    
    # remove tick labels 
    ax1.set_yticklabels([])
    ax1.set_xticklabels([])
    ax1.set_zticklabels([])
        
    # legend
    ax1.legend(loc=(0.2, -0.08), fancybox=True, facecolor='white', edgecolor='black', frameon=True)
    
    # set title
    ax1.set_title('Three body orbits (Classical Solution)', fontsize=14, color = 'white')
        
    # plot x_2 post newtonian

    ax2 = fig.add_subplot(1, 2, 2, projection='3d')
    #ax1.figure.clear()
    ax2.grid(False)

    line4, = ax2.plot([], [], color = '#A325D9', linewidth=0.8)
    line5, = ax2.plot([], [], color = '#04AD7B', linewidth=0.8)
    line6, = ax2.plot([], [], color = '#D93425', linewidth=0.8)
    
    point4, = ax2.plot([], [], marker='o', color = '#A325D9', markersize=4)
    point5, = ax2.plot([], [], marker='o', color = '#04AD7B', markersize=4)
    point6, = ax2.plot([], [], marker='o', color = '#D93425', markersize=4)

    ax2.set_facecolor('#021226')
    ax2.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.2))
    ax2.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.2))
    ax2.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.2))
    # Take out the axis
    #ax1.set_axis_off()
    
    # set point-of-view: specified by (altitude degrees, azimuth degree)
    ax2.view_init(elev=23, azim=-131)
    
    # set limits 
    ax2.set_xlim(-2, 2)
    ax2.set_ylim(-2, 2)
    ax2.set_zlim(-2, 2)
    
    # # set axis labels
    # ax1.set_xlabel('X', fontsize=10)
    # ax1.set_ylabel('Y', fontsize=10)
    # ax1.set_zlabel('Z', fontsize=10)
    
    # remove tick labels 
    ax2.set_yticklabels([])
    ax2.set_xticklabels([])
    ax2.set_zticklabels([])
        
    # legend
    #ax2.text(0.2, -0.08, s= "Feel free to rotate the plot!")
    
    # set title
    ax2.set_title('Three body orbits (Post-Newtonian solution)', fontsize=14, color = 'white')


    def update(i):
        i = 250*i
        if i >= 150000:
            i = 150000-1
        xq_1 = x_1[0:i]
        yq_1 = y_1[0:i]
        zq_1 = z_1[0:i]

        xq_2 = x_2[0:i]
        yq_2 = y_2[0:i]
        zq_2 = z_2[0:i]
        
        xq_3 = x_3[0:i]
        yq_3 = y_3[0:i]
        zq_3 = z_3[0:i]

        xq_4 = x_4[0:i]
        yq_4 = y_4[0:i]
        zq_4 = z_4[0:i]

        xq_5 = x_5[0:i]
        yq_5 = y_5[0:i]
        zq_5 = z_5[0:i]
        
        xq_6 = x_6[0:i]
        yq_6 = y_6[0:i]
        zq_6 = z_6[0:i]

        tq = t[0:i]
        
        line.set_data(xq_1, yq_1)
        line.set_3d_properties(zq_1)
        point.set_data(x_1[i], y_1[i])
        point.set_3d_properties(z_1[i])
        
        line2.set_data(xq_2, yq_2)
        line2.set_3d_properties(zq_2)
        point2.set_data(x_2[i], y_2[i])
        point2.set_3d_properties(z_2[i])
        
        line3.set_data(xq_3, yq_3)
        line3.set_3d_properties(zq_3)
        point3.set_data(x_3[i], y_3[i])
        point3.set_3d_properties(z_3[i])

        line4.set_data(xq_4, yq_4)
        line4.set_3d_properties(zq_4)
        point4.set_data(x_4[i], y_4[i])
        point4.set_3d_properties(z_4[i])
        
        line5.set_data(xq_5, yq_5)
        line5.set_3d_properties(zq_5)
        point5.set_data(x_5[i], y_5[i])
        point5.set_3d_properties(z_5[i])
        
        line6.set_data(xq_6, yq_6)
        line6.set_3d_properties(zq_6)
        point6.set_data(x_6[i], y_6[i])
        point6.set_3d_properties(z_6[i])
        
        return line, line2, line3, line4, line5, line6, point, point2, point3, point4, point5, point6
    
    # instantiate the animator
    global ani
    ani = FuncAnimation(fig, update, frames=np.size(x_1), interval=0, blit=True)
    #ani.event_source.stop()
    # save as mp4. This requires mplayer or ffmpeg to be installed
    # anim.save('RosslerAttractor.mp4', fps=15, extra_args=['-vcodec', 'libx264'])
    
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
        button3 = QtWidgets.QPushButton("Figure 8! (m1 = m2 = m3 = 1, so the mass slider does not affects the system)")
        button3.setGeometry(10, 7, 10, 3)
        button3.clicked.connect(self.figure_8)
        button4 = QtWidgets.QPushButton("Circle! (m1 = m2 = m3, so only the mass slider 1 affects the system)")
        button4.clicked.connect(self.circle)
        combobox = QtWidgets.QComboBox(self)
        # Add preconfigured cases to the box
        combobox.addItem("Figure 8")
        combobox.addItem("Circle")
        combobox.addItem("M0 (Triquette)")
        combobox.addItem("Planar M1 (Bow Tie)")
        combobox.addItem("Planar M2 (Butterfly)")
        combobox.addItem("Planar M3")
        combobox.addItem("Planar Unstable1")
        combobox.addItem("Random Configuration")
        choose = QtWidgets.QLabel("Choose a predefined configuration or have fun with random orbits! (Don't forget to stop the simulation before trying another config).")
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
        groupBox.setTitle('McGill Physics Hackathon 2021: Three-Body Gravitational Problem')
        groupBox.setStyleSheet(self.getStyleSheet("./styles.css"))
        groupBox.setAlignment(QtCore.Qt.AlignCenter)
        lay.addLayout(vlay1)
        lay.addLayout(vlay2)
        lay.addLayout(vlay3)
        lay.addLayout(vlay4)
        lay.addLayout(vlay5)
        lay.addLayout(vlay6)
        groupBox.setLayout(lay)

        choose2 = QtWidgets.QLabel("You can also choose the masses and run a custom built simulation with random positions and velocities.")
        
        # Sliders for the masses
        a_label = QtWidgets.QLabel("Massa João =")
        a_slider = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        a_slider.setRange(0,50)
        a_slider.setValue(1)
        self.a_slider = a_slider

        button = QtWidgets.QPushButton("Ready, Set, Go!")
        button.clicked.connect(self.main)

        vlay3.addWidget(choose2)
        vlay5.addWidget(button)
        vlay5.addWidget(button2)
        b_label = QtWidgets.QLabel("Massa Olga = ")
        b_slider = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        b_slider.setRange(0,50)
        b_slider.setValue(1)
        self.b_slider = b_slider

        c_label = QtWidgets.QLabel("Massa Pedrin = ")
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
        self.spinBoxm1 = QtWidgets.QDoubleSpinBox()
        self.spinBoxm1.setRange(0, 50)
        self.spinBoxm1.setValue(1)
        self.spinBoxm1.setPrefix("m1 =  ")
        self.spinBoxm1.valueChanged.connect(a_value.setNum)
        vlay4.addWidget(self.spinBoxm1)

        self.spinBoxm2 = QtWidgets.QDoubleSpinBox()
        self.spinBoxm2.setRange(0, 50)
        self.spinBoxm2.setValue(1)
        self.spinBoxm2.setPrefix("m2 =  ")
        self.spinBoxm2.valueChanged.connect(b_value.setNum)
        vlay4.addWidget(self.spinBoxm2)

        self.spinBoxm3 = QtWidgets.QDoubleSpinBox()
        self.spinBoxm3.setRange(0, 50)
        self.spinBoxm3.setValue(2)
        self.spinBoxm3.setPrefix("m3 =  ")
        self.spinBoxm3.valueChanged.connect(c_value.setNum)
        vlay4.addWidget(self.spinBoxm3)

        #lay.addWidget(button2)
        #self.plot()


    def case(self, text):
        if text == "Figure 8":
            self.figure_8()
        elif text == "Circle":
            self.circle()   
        elif text == "M0 (Triquette)":
            self.triquette()   
        elif text == "Planar Unstable1":
            self.U1()          
        elif text == "Planar M1 (Bow Tie)":
            self.M1()   
        elif text == "Planar M2 (Butterfly)":
            self.M2() 
        elif text == "Planar M3":
            self.M3() 
        elif text == "Random Configuration":
            self.main()

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

    def U1(self):
        #ani.event_source.stop()
        self.figure.clear()
        m1 = 1
        m2 = m1
        m3 = 2
        # time interval and the step size
        t = np.arange(0, 30, anispeed)
        
        # vectors for the solutions
        x_1 = np.zeros((len(t)))
        y_1 = np.zeros((len(t)))
        z_1 = np.zeros((len(t)))
        
        x_2 = np.zeros((len(t)))
        y_2 = np.zeros((len(t)))
        z_2 = np.zeros((len(t)))

        x_3 = np.zeros((len(t)))
        y_3 = np.zeros((len(t)))
        z_3 = np.zeros((len(t)))       

        x_4 = np.zeros((len(t)))
        y_4 = np.zeros((len(t)))
        z_4 = np.zeros((len(t)))
        
        x_5 = np.zeros((len(t)))
        y_5 = np.zeros((len(t)))
        z_5 = np.zeros((len(t)))

        x_6 = np.zeros((len(t)))
        y_6 = np.zeros((len(t)))
        z_6 = np.zeros((len(t)))

        # initial conditions
        """
        This initial conditions generate the triquette figure
        """
        L = 2
        #Define initial position vectors
        r1=[-1,0,0.0]
        r2=[1,0,0.0]
        r3=[0.001,0,0.0]
        Pos_0 = np.array([r1,r2,r3])
        #Convert pos vectors to arrays
        r1=np.array(r1,dtype="float64")
        r2=np.array(r2,dtype="float64")
        r3=np.array(r3,dtype="float64")

        #Define initial velocities
        v1=[0.305722433, 0.521512426, 0.0]
        v2=v1
        v3=[-2*v1[0]*m1/m3, -2*v1[1]*m2/m3, 0.0]
        Vel_0 = np.array([v1,v2,v3])
        #Convert velocity vectors to arrays
        v1=np.array(v1,dtype="float64")
        v2=np.array(v2,dtype="float64")
        v3=np.array(v3,dtype="float64")
        #Update COM formula
        r_com=(m1*r1+m2*r2+m3*r3)/(m1+m2+m3)
        #Update velocity of COM formula
        v_com=(m1*v1+m2*v2+m3*v3)/(m1+m2+m3)
        
        init_params=np.array([r1,r2,r3,v1,v2,v3]) #Initial parameters
        init_params=init_params.flatten() #Flatten to make 1D array
        
        #three_body_sol=sp.integrate.odeint(ThreeBodyEquations,init_params,t,args=(G,m1,m2,m3))
        
        three_body_sol=solve_ivp(ThreeBodyEquations, [0, 30], init_params, t_eval = t, method = 'Radau', max_step = 1e-2, args=(G,m1,m2,m3))

        r1_sol=three_body_sol.y[:3]
        r2_sol=three_body_sol.y[3:6]
        r3_sol=three_body_sol.y[6:9]
        
        x_1 = r1_sol[0]
        y_1 = r1_sol[1]
        z_1 = r1_sol[2]

        x_2 = r2_sol[0]
        y_2 = r2_sol[1]
        z_2 = r2_sol[2]
        
        x_3 = r3_sol[0]
        y_3 = r3_sol[1]
        z_3 = r3_sol[2]
        
        # Post Newtonian Solution

        u0 = np.zeros(18)

        u0[:3] = Pos_0[0,:]
        u0[3:6] = Vel_0[0,:]
        u0[6:9] = Pos_0[1,:]
        u0[9:12] = Vel_0[1,:]
        u0[12:15] = Pos_0[2,:]
        u0[15:] = Vel_0[2, :]
        
        three_body_sol=solve_ivp(dydx, [0, 30], u0, t_eval = t, method = 'Radau', max_step = 1e-2, args=(m1,m2,m3))

        r4_sol=three_body_sol.y[:3]
        r5_sol=three_body_sol.y[6:9]
        r6_sol=three_body_sol.y[12:15]
        
        x_4 = r4_sol[0]
        y_4 = r4_sol[1]
        z_4 = r4_sol[2]

        x_5 = r5_sol[0]
        y_5 = r5_sol[1]
        z_5 = r5_sol[2]
        
        x_6 = r6_sol[0]
        y_6 = r6_sol[1]
        z_6 = r6_sol[2]
        
        plot_solution(self, fig, x_1, y_1, z_1, t, x_2, y_2, z_2, x_3, y_3, z_3, 
         x_4, y_4, z_4, x_5, y_5, z_5, x_6, y_6, z_6)
        
    def M3(self):
        #ani.event_source.stop()
        self.figure.clear()
        m1 = 1
        m2 = m1
        m3 = 0.5
        # time interval and the step size
        t = np.arange(0, 30, anispeed)
        
        # vectors for the solutions
        x_1 = np.zeros((len(t)))
        y_1 = np.zeros((len(t)))
        z_1 = np.zeros((len(t)))
        
        x_2 = np.zeros((len(t)))
        y_2 = np.zeros((len(t)))
        z_2 = np.zeros((len(t)))

        x_3 = np.zeros((len(t)))
        y_3 = np.zeros((len(t)))
        z_3 = np.zeros((len(t)))       

        x_4 = np.zeros((len(t)))
        y_4 = np.zeros((len(t)))
        z_4 = np.zeros((len(t)))
        
        x_5 = np.zeros((len(t)))
        y_5 = np.zeros((len(t)))
        z_5 = np.zeros((len(t)))

        x_6 = np.zeros((len(t)))
        y_6 = np.zeros((len(t)))
        z_6 = np.zeros((len(t)))

        # initial conditions
        """
        This initial conditions generate the triquette figure
        """
        L = 2
        #Define initial position vectors
        r1=[-1,0,0.0]
        r2=[1,0,0.0]
        r3=[0.00001,0,0.0]
        Pos_0=np.array([r1,r2,r3])
        #Convert pos vectors to arrays
        r1=np.array(r1,dtype="float64")
        r2=np.array(r2,dtype="float64")
        r3=np.array(r3,dtype="float64")
        #Define initial velocities
        v1=[0.3420307307, 0.1809369236, 0.0]
        v2=v1
        v3=[-2*v1[0]*m1/m3, -2*v1[1]*m2/m3, 0.0]
        Vel_0=np.array([v1,v2,v3])
        #Convert velocity vectors to arrays
        v1=np.array(v1,dtype="float64")
        v2=np.array(v2,dtype="float64")
        v3=np.array(v3,dtype="float64")
        
        #Update COM formula
        r_com=(m1*r1+m2*r2+m3*r3)/(m1+m2+m3)
        #Update velocity of COM formula
        v_com=(m1*v1+m2*v2+m3*v3)/(m1+m2+m3)
        
        init_params=np.array([r1,r2,r3,v1,v2,v3]) #Initial parameters
        init_params=init_params.flatten() #Flatten to make 1D array
        
        #three_body_sol=sp.integrate.odeint(ThreeBodyEquations,init_params,t,args=(G,m1,m2,m3))
        
        three_body_sol=solve_ivp(ThreeBodyEquations, [0, 30], init_params, t_eval = t, method = 'Radau', max_step = 1e-2, args=(G,m1,m2,m3))

        r1_sol=three_body_sol.y[:3]
        r2_sol=three_body_sol.y[3:6]
        r3_sol=three_body_sol.y[6:9]
        
        x_1 = r1_sol[0]
        y_1 = r1_sol[1]
        z_1 = r1_sol[2]

        x_2 = r2_sol[0]
        y_2 = r2_sol[1]
        z_2 = r2_sol[2]
        
        x_3 = r3_sol[0]
        y_3 = r3_sol[1]
        z_3 = r3_sol[2]
        
        # Post Newtonian Solution
        u0 = np.zeros(18)

        u0[:3] = Pos_0[0,:]
        u0[3:6] = Vel_0[0,:]
        u0[6:9] = Pos_0[1,:]
        u0[9:12] = Vel_0[1,:]
        u0[12:15] = Pos_0[2,:]
        u0[15:] = Vel_0[2, :]
        
        three_body_sol=solve_ivp(dydx, [0, 30], u0, t_eval = t, method = 'Radau', max_step = 1e-2, args=(m1,m2,m3))

        r4_sol=three_body_sol.y[:3]
        r5_sol=three_body_sol.y[6:9]
        r6_sol=three_body_sol.y[12:15]
        
        x_4 = r4_sol[0]
        y_4 = r4_sol[1]
        z_4 = r4_sol[2]

        x_5 = r5_sol[0]
        y_5 = r5_sol[1]
        z_5 = r5_sol[2]
        
        x_6 = r6_sol[0]
        y_6 = r6_sol[1]
        z_6 = r6_sol[2]
        
        plot_solution(self, fig, x_1, y_1, z_1, t, x_2, y_2, z_2, x_3, y_3, z_3, 
         x_4, y_4, z_4, x_5, y_5, z_5, x_6, y_6, z_6)
        
    def M2(self):
        #ani.event_source.stop()
        self.figure.clear()
        m1 = 1
        m2 = m1
        m3 = 2
        # time interval and the step size
        t = np.arange(0, 30, anispeed)

        # vectors for the solutions
        x_1 = np.zeros((len(t)))
        y_1 = np.zeros((len(t)))
        z_1 = np.zeros((len(t)))
        
        x_2 = np.zeros((len(t)))
        y_2 = np.zeros((len(t)))
        z_2 = np.zeros((len(t)))

        x_3 = np.zeros((len(t)))
        y_3 = np.zeros((len(t)))
        z_3 = np.zeros((len(t)))       

        x_4 = np.zeros((len(t)))
        y_4 = np.zeros((len(t)))
        z_4 = np.zeros((len(t)))
        
        x_5 = np.zeros((len(t)))
        y_5 = np.zeros((len(t)))
        z_5 = np.zeros((len(t)))

        x_6 = np.zeros((len(t)))
        y_6 = np.zeros((len(t)))
        z_6 = np.zeros((len(t)))

        
        # initial conditions
        """
        This initial conditions generate the triquette figure
        """
        L = 2
        #Define initial position vectors
        r1=[-1,0,0.0]
        r2=[1,0,0.0]
        r3=[0.00001,0,0.0]
        Pos_0=np.array([r1,r2,r3])
        #Convert pos vectors to arrays
        r1=np.array(r1,dtype="float64")
        r2=np.array(r2,dtype="float64")
        r3=np.array(r3,dtype="float64")
        #Define initial velocities
        v1=[0.664910758, 0.832416786, 0.0]
        v2=v1
        v3=[-2*v1[0]*m1/m3, -2*v1[1]*m2/m3, 0.0]
        Vel_0=np.array([v1,v2,v3])
        #Convert velocity vectors to arrays
        v1=np.array(v1,dtype="float64")
        v2=np.array(v2,dtype="float64")
        v3=np.array(v3,dtype="float64")
        
        #Update COM formula
        r_com=(m1*r1+m2*r2+m3*r3)/(m1+m2+m3)
        #Update velocity of COM formula
        v_com=(m1*v1+m2*v2+m3*v3)/(m1+m2+m3)
        
        init_params=np.array([r1,r2,r3,v1,v2,v3]) #Initial parameters
        init_params=init_params.flatten() #Flatten to make 1D array
        
        #three_body_sol=sp.integrate.odeint(ThreeBodyEquations,init_params,t,args=(G,m1,m2,m3))
        
        three_body_sol=solve_ivp(ThreeBodyEquations, [0, 30], init_params, t_eval = t, method = 'Radau', max_step = 1e-2, args=(G,m1,m2,m3))

        r1_sol=three_body_sol.y[:3]
        r2_sol=three_body_sol.y[3:6]
        r3_sol=three_body_sol.y[6:9]
        
        x_1 = r1_sol[0]
        y_1 = r1_sol[1]
        z_1 = r1_sol[2]

        x_2 = r2_sol[0]
        y_2 = r2_sol[1]
        z_2 = r2_sol[2]
        
        x_3 = r3_sol[0]
        y_3 = r3_sol[1]
        z_3 = r3_sol[2]
        
        # Post Newtonian Solution
        u0 = np.zeros(18)

        u0[:3] = Pos_0[0,:]
        u0[3:6] = Vel_0[0,:]
        u0[6:9] = Pos_0[1,:]
        u0[9:12] = Vel_0[1,:]
        u0[12:15] = Pos_0[2,:]
        u0[15:] = Vel_0[2, :]
        
        three_body_sol=solve_ivp(dydx, [0, 30], u0, t_eval = t, method = 'Radau', max_step = 1e-2, args=(m1,m2,m3))

        r4_sol=three_body_sol.y[:3]
        r5_sol=three_body_sol.y[6:9]
        r6_sol=three_body_sol.y[12:15]
        
        x_4 = r4_sol[0]
        y_4 = r4_sol[1]
        z_4 = r4_sol[2]

        x_5 = r5_sol[0]
        y_5 = r5_sol[1]
        z_5 = r5_sol[2]
        
        x_6 = r6_sol[0]
        y_6 = r6_sol[1]
        z_6 = r6_sol[2]
        
        plot_solution(self, fig, x_1, y_1, z_1, t, x_2, y_2, z_2, x_3, y_3, z_3, 
         x_4, y_4, z_4, x_5, y_5, z_5, x_6, y_6, z_6)
        
    def M1(self):
        #ani.event_source.stop()
        self.figure.clear()
        m1 = 1
        m2 = m1
        m3 = 4
        # time interval and the step size
        t = np.arange(0, 30, anispeed)
        
        # vectors for the solutions
        x_1 = np.zeros((len(t)))
        y_1 = np.zeros((len(t)))
        z_1 = np.zeros((len(t)))
        
        x_2 = np.zeros((len(t)))
        y_2 = np.zeros((len(t)))
        z_2 = np.zeros((len(t)))

        x_3 = np.zeros((len(t)))
        y_3 = np.zeros((len(t)))
        z_3 = np.zeros((len(t)))       

        x_4 = np.zeros((len(t)))
        y_4 = np.zeros((len(t)))
        z_4 = np.zeros((len(t)))
        
        x_5 = np.zeros((len(t)))
        y_5 = np.zeros((len(t)))
        z_5 = np.zeros((len(t)))

        x_6 = np.zeros((len(t)))
        y_6 = np.zeros((len(t)))
        z_6 = np.zeros((len(t)))

        # initial conditions
        """
        This initial conditions generate the triquette figure
        """
        L = 2
        #Define initial position vectors
        r1=[-1,0.0,0.0]
        r2=[1.0,0.0,0.0]
        r3=[0.0001,0.0,0.0]
        Pos_0=np.array([r1,r2,r3])
        #Convert pos vectors to arrays
        r1=np.array(r1,dtype="float64")
        r2=np.array(r2,dtype="float64")
        r3=np.array(r3,dtype="float64")
        #Define initial velocities
        v1=[0.991198122, 0.711947212, 0.0]
        v2=v1
        v3=[-2*v1[0]*m1/m3, -2*v1[1]*m2/m3, 0.0]
        Vel_0=np.array([v1,v2,v3])
        #Convert velocity vectors to arrays
        v1=np.array(v1,dtype="float64")
        v2=np.array(v2,dtype="float64")
        v3=np.array(v3,dtype="float64")

        #Update COM formula
        r_com=(m1*r1+m2*r2+m3*r3)/(m1+m2+m3)
        #Update velocity of COM formula
        v_com=(m1*v1+m2*v2+m3*v3)/(m1+m2+m3)
        
        init_params=np.array([r1,r2,r3,v1,v2,v3]) #Initial parameters
        init_params=init_params.flatten() #Flatten to make 1D array
        
        #three_body_sol=sp.integrate.odeint(ThreeBodyEquations,init_params,t,args=(G,m1,m2,m3))
        
        three_body_sol=solve_ivp(ThreeBodyEquations, [0, 30], init_params, t_eval = t, method = 'Radau', max_step = 1e-2, args=(G,m1,m2,m3))

        r1_sol=three_body_sol.y[:3]
        r2_sol=three_body_sol.y[3:6]
        r3_sol=three_body_sol.y[6:9]
        
        x_1 = r1_sol[0]
        y_1 = r1_sol[1]
        z_1 = r1_sol[2]

        x_2 = r2_sol[0]
        y_2 = r2_sol[1]
        z_2 = r2_sol[2]
        
        x_3 = r3_sol[0]
        y_3 = r3_sol[1]
        z_3 = r3_sol[2]
        
        # Post Newtonian Solution
        u0 = np.zeros(18)

        u0[:3] = Pos_0[0,:]
        u0[3:6] = Vel_0[0,:]
        u0[6:9] = Pos_0[1,:]
        u0[9:12] = Vel_0[1,:]
        u0[12:15] = Pos_0[2,:]
        u0[15:] = Vel_0[2, :]
        
        three_body_sol=solve_ivp(dydx, [0, 30], u0, t_eval = t, method = 'Radau', max_step = 1e-2, args=(m1,m2,m3))

        r4_sol=three_body_sol.y[:3]
        r5_sol=three_body_sol.y[6:9]
        r6_sol=three_body_sol.y[12:15]
        
        x_4 = r4_sol[0]
        y_4 = r4_sol[1]
        z_4 = r4_sol[2]

        x_5 = r5_sol[0]
        y_5 = r5_sol[1]
        z_5 = r5_sol[2]
        
        x_6 = r6_sol[0]
        y_6 = r6_sol[1]
        z_6 = r6_sol[2]
        
        plot_solution(self, fig, x_1, y_1, z_1, t, x_2, y_2, z_2, x_3, y_3, z_3, 
         x_4, y_4, z_4, x_5, y_5, z_5, x_6, y_6, z_6)
        
    def triquette(self):
        #ani.event_source.stop()
        self.figure.clear()
        m1 = 2.5
        m2 = m1
        m3 = m1
        # time interval and the step size
        t = np.arange(0, 30, anispeed)
        
        # vectors for the solutions
        x_1 = np.zeros((len(t)))
        y_1 = np.zeros((len(t)))
        z_1 = np.zeros((len(t)))
        
        x_2 = np.zeros((len(t)))
        y_2 = np.zeros((len(t)))
        z_2 = np.zeros((len(t)))
        
        # initial conditions
        """
        This initial conditions generate the triquette figure
        """
        L = 5
        #Define initial position vectors
        r1=[-L/2,-((np.sqrt(3)*L/2) - L/np.sqrt(3)),0.0]
        r2=[0.0,(np.sqrt(3)*L/2)-((np.sqrt(3)*L/2) - L/np.sqrt(3)),0.0]
        r3=[L/2,-((np.sqrt(3)*L/2) -L/np.sqrt(3)),0.0]

        #Convert pos vectors to arrays
        r1=np.array(r1,dtype="float64")
        r2=np.array(r2,dtype="float64")
        r3=np.array(r3,dtype="float64")
        Pos_0 = np.array([r1,r2,r3])
        #Define initial velocities
        v1=[np.sqrt(G/L)*1/2, -np.sqrt(G/L)*np.sqrt(3)/2, 0.0]
        v2=[-np.sqrt(G/L), 0.0,0.0]
        v3=[np.sqrt(G/L)*1/2, np.sqrt(G/L)*np.sqrt(3)/2, 0.0]
        Vel_0 = np.array([v1,v2,v3])
        #Convert velocity vectors to arrays
        v1=np.array(v1,dtype="float64")
        v2=np.array(v2,dtype="float64")
        v3=np.array(v3,dtype="float64")
        
        # #Update COM formula
        # r_com=(m1*r1+m2*r2+m3*r3)/(m1+m2+m3)
        # #Update velocity of COM formula
        # v_com=(m1*v1+m2*v2+m3*v3)/(m1+m2+m3)
        
        init_params=np.array([r1,r2,r3,v1,v2,v3]) #Initial parameters
        init_params=init_params.flatten() #Flatten to make 1D array
        
        #three_body_sol=sp.integrate.odeint(ThreeBodyEquations,init_params,t,args=(G,m1,m2,m3))
        
        three_body_sol=solve_ivp(ThreeBodyEquations, [0, 30], init_params, t_eval = t, method = 'Radau', max_step = 1e-2, args=(G,m1,m2,m3))

        r1_sol=three_body_sol.y[:3]
        r2_sol=three_body_sol.y[3:6]
        r3_sol=three_body_sol.y[6:9]
        
        x_1 = r1_sol[0]
        y_1 = r1_sol[1]
        z_1 = r1_sol[2]

        x_2 = r2_sol[0]
        y_2 = r2_sol[1]
        z_2 = r2_sol[2]
        
        x_3 = r3_sol[0]
        y_3 = r3_sol[1]
        z_3 = r3_sol[2]

        x_4 = np.zeros((len(t)))
        y_4 = np.zeros((len(t)))
        z_4 = np.zeros((len(t)))
        
        x_5 = np.zeros((len(t)))
        y_5 = np.zeros((len(t)))
        z_5 = np.zeros((len(t)))

        x_6 = np.zeros((len(t)))
        y_6 = np.zeros((len(t)))
        z_6 = np.zeros((len(t)))

        # Post Newtonian Solution
        u0 = np.zeros(18)

        u0[:3] = Pos_0[0,:]
        u0[3:6] = Vel_0[0,:]
        u0[6:9] = Pos_0[1,:]
        u0[9:12] = Vel_0[1,:]
        u0[12:15] = Pos_0[2,:]
        u0[15:] = Vel_0[2, :]
        
        three_body_sol=solve_ivp(dydx, [0, 30], u0, t_eval = t, method = 'Radau', max_step = 1e-2, args=(m1,m2,m3))

        r4_sol=three_body_sol.y[:3]
        r5_sol=three_body_sol.y[6:9]
        r6_sol=three_body_sol.y[12:15]
        
        x_4 = r4_sol[0]
        y_4 = r4_sol[1]
        z_4 = r4_sol[2]

        x_5 = r5_sol[0]
        y_5 = r5_sol[1]
        z_5 = r5_sol[2]
        
        x_6 = r6_sol[0]
        y_6 = r6_sol[1]
        z_6 = r6_sol[2]
        
        plot_solution(self, fig, x_1, y_1, z_1, t, x_2, y_2, z_2, x_3, y_3, z_3, 
         x_4, y_4, z_4, x_5, y_5, z_5, x_6, y_6, z_6)
        
    def figure_8(self):
        # initial masses
        m1 = 1
        m2 = 1
        m3 = 1
        # time interval and the step size
        t = np.arange(0, 30, anispeed)
        # vectors for the solutions
        x_1 = np.zeros((len(t)))
        y_1 = np.zeros((len(t)))
        z_1 = np.zeros((len(t)))
        
        x_2 = np.zeros((len(t)))
        y_2 = np.zeros((len(t)))
        z_2 = np.zeros((len(t)))

        x_3 = np.zeros((len(t)))
        y_3 = np.zeros((len(t)))
        z_3 = np.zeros((len(t)))       

        x_4 = np.zeros((len(t)))
        y_4 = np.zeros((len(t)))
        z_4 = np.zeros((len(t)))
        
        x_5 = np.zeros((len(t)))
        y_5 = np.zeros((len(t)))
        z_5 = np.zeros((len(t)))

        x_6 = np.zeros((len(t)))
        y_6 = np.zeros((len(t)))
        z_6 = np.zeros((len(t)))

        # initial conditions
        """
        This initial conditions generate the eight figure
        """
        #Define initial position vectors
        r1=[-0.97000436, 0.24308753,0]
        r2=[0,0,0]
        r3=[0.97000436, -0.24308753,0]
        #Convert pos vectors to arrays
        r1=np.array(r1,dtype="float64")
        r2=np.array(r2,dtype="float64")
        r3=np.array(r3,dtype="float64")
        Pos_0 = np.array([r1,r2,r3])
        #Define initial velocities
        v1=[0.4662036850, 0.4323657300,0]
        v2=[-0.93240737, -0.86473146,0]
        v3=[0.4662036850, 0.4323657300,0]
        Vel_0 = np.array([v1,v2,v3])
        #Convert velocity vectors to arrays
        v1=np.array(v1,dtype="float64")
        v2=np.array(v2,dtype="float64")
        v3=np.array(v3,dtype="float64")
        # #Update COM formula
        # r_com=(m1*r1+m2*r2+m3*r3)/(m1+m2+m3)
        # #Update velocity of COM formula
        # v_com=(m1*v1+m2*v2+m3*v3)/(m1+m2+m3)
        
        init_params=np.array([r1,r2,r3,v1,v2,v3]) #Initial parameters
        init_params=init_params.flatten() #Flatten to make 1D array
        
        three_body_sol=solve_ivp(ThreeBodyEquations, [0, 30], init_params, t_eval = t, method = 'Radau', max_step = 1e-2, args=(G,m1,m2,m3))
        #three_body_sol=sp.integrate.odeint(ThreeBodyEquations,init_params,t,args=(G,m1,m2,m3))
        
        r1_sol=three_body_sol.y[:3]
        r2_sol=three_body_sol.y[3:6]
        r3_sol=three_body_sol.y[6:9]
        
        x_1 = r1_sol[0]
        y_1 = r1_sol[1]
        z_1 = r1_sol[2]

        x_2 = r2_sol[0]
        y_2 = r2_sol[1]
        z_2 = r2_sol[2]
        
        x_3 = r3_sol[0]
        y_3 = r3_sol[1]
        z_3 = r3_sol[2]
        
        # Post Newtonian Solution

        u0 = np.zeros(18)

        u0[:3] = Pos_0[0,:]
        u0[3:6] = Vel_0[0,:]
        u0[6:9] = Pos_0[1,:]
        u0[9:12] = Vel_0[1,:]
        u0[12:15] = Pos_0[2,:]
        u0[15:] = Vel_0[2, :]
        
        three_body_sol=solve_ivp(dydx, [0, 30], u0, t_eval = t, method = 'Radau', max_step = 1e-2, args=(m1,m2,m3))

        r4_sol=three_body_sol.y[:3]
        r5_sol=three_body_sol.y[6:9]
        r6_sol=three_body_sol.y[12:15]
        
        x_4 = r4_sol[0]
        y_4 = r4_sol[1]
        z_4 = r4_sol[2]

        x_5 = r5_sol[0]
        y_5 = r5_sol[1]
        z_5 = r5_sol[2]
        
        x_6 = r6_sol[0]
        y_6 = r6_sol[1]
        z_6 = r6_sol[2]
        
        plot_solution(self, fig, x_1, y_1, z_1, t, x_2, y_2, z_2, x_3, y_3, z_3, 
         x_4, y_4, z_4, x_5, y_5, z_5, x_6, y_6, z_6)
        
    def circle(self):
        #ani.event_source.stop()
        self.figure.clear()
        m1 = self.a_slider.value()
        m2 = m1
        m3 = m1
        # time interval and the step size
        t = np.arange(0, 30, anispeed)
        
        # vectors for the solutions
        x_1 = np.zeros((len(t)))
        y_1 = np.zeros((len(t)))
        z_1 = np.zeros((len(t)))
        
        x_2 = np.zeros((len(t)))
        y_2 = np.zeros((len(t)))
        z_2 = np.zeros((len(t)))
        
        x_4 = np.zeros((len(t)))
        y_4 = np.zeros((len(t)))
        z_4 = np.zeros((len(t)))
        
        x_5 = np.zeros((len(t)))
        y_5 = np.zeros((len(t)))
        z_5 = np.zeros((len(t)))

        x_6 = np.zeros((len(t)))
        y_6 = np.zeros((len(t)))
        z_6 = np.zeros((len(t)))

        # initial conditions
        """
        This initial conditions generate the circle
        """
        L = 5
        #Define initial position vectors
        r1=[-L/2,-((np.sqrt(3)*L/2) - L/np.sqrt(3)),0.0]
        r2=[0.0,(np.sqrt(3)*L/2)-((np.sqrt(3)*L/2) - L/np.sqrt(3)),0.0]
        r3=[L/2,-((np.sqrt(3)*L/2) -L/np.sqrt(3)),0.0]
        Pos_0 = np.array([r1,r2,r3])
        #Convert pos vectors to arrays
        r1=np.array(r1,dtype="float64")
        r2=np.array(r2,dtype="float64")
        r3=np.array(r3,dtype="float64")
        #Define initial velocities
        v1=[np.sqrt(G*m1/L)*1/2, -np.sqrt(G*m1/L)*np.sqrt(3)/2, 0.0]
        v2=[-1*np.sqrt(G*m2/L),0.0*np.sqrt(G*m2/L),0.0]
        v3=[np.sqrt(G*m3/L)*1/2, np.sqrt(G*m3/L)*np.sqrt(3)/2, 0.0]
        Vel_0 = np.array([v1,v2,v3])
        #Convert velocity vectors to arrays
        v1=np.array(v1,dtype="float64")
        v2=np.array(v2,dtype="float64")
        v3=np.array(v3,dtype="float64")
        
        #Update COM formula
        r_com=(m1*r1+m2*r2+m3*r3)/(m1+m2+m3)
        #Update velocity of COM formula
        v_com=(m1*v1+m2*v2+m3*v3)/(m1+m2+m3)
        
        init_params=np.array([r1,r2,r3,v1,v2,v3]) #Initial parameters
        init_params=init_params.flatten() #Flatten to make 1D array
        
        three_body_sol=solve_ivp(ThreeBodyEquations, [0, 30], init_params, t_eval = t, method = 'Radau', max_step = 1e-2, args=(G,m1,m2,m3))
        #three_body_sol=sp.integrate.odeint(ThreeBodyEquations,init_params,t,args=(G,m1,m2,m3))
        
        r1_sol=three_body_sol.y[:3]
        r2_sol=three_body_sol.y[3:6]
        r3_sol=three_body_sol.y[6:9]
        
        x_1 = r1_sol[0]
        y_1 = r1_sol[1]
        z_1 = r1_sol[2]

        x_2 = r2_sol[0]
        y_2 = r2_sol[1]
        z_2 = r2_sol[2]
        
        x_3 = r3_sol[0]
        y_3 = r3_sol[1]
        z_3 = r3_sol[2]

       # Post Newtonian Solution

        u0 = np.zeros(18)

        u0[:3] = Pos_0[0,:]
        u0[3:6] = Vel_0[0,:]
        u0[6:9] = Pos_0[1,:]
        u0[9:12] = Vel_0[1,:]
        u0[12:15] = Pos_0[2,:]
        u0[15:] = Vel_0[2, :]
        
        three_body_sol=solve_ivp(dydx, [0, 30], u0, t_eval = t, method = 'Radau', max_step = 1e-2, args=(m1,m2,m3))

        r4_sol=three_body_sol.y[:3]
        r5_sol=three_body_sol.y[6:9]
        r6_sol=three_body_sol.y[12:15]
        
        x_4 = r4_sol[0]
        y_4 = r4_sol[1]
        z_4 = r4_sol[2]

        x_5 = r5_sol[0]
        y_5 = r5_sol[1]
        z_5 = r5_sol[2]
        
        x_6 = r6_sol[0]
        y_6 = r6_sol[1]
        z_6 = r6_sol[2]
        
        plot_solution(self, fig, x_1, y_1, z_1, t, x_2, y_2, z_2, x_3, y_3, z_3, 
         x_4, y_4, z_4, x_5, y_5, z_5, x_6, y_6, z_6)
           
        
        #plot_solution(self, fig, x_1, y_1, z_1, t, x_2, y_2, z_2, x_3, y_3, z_3)
        
    def main(self):
        #ani.event_source.stop()
        self.figure.clear()
        m1 = self.spinBoxm1.value()  
        m2 = self.spinBoxm2.value()
        m3 = self.spinBoxm3.value()
        # time interval and the step size
        t = np.arange(0, 30, anispeed)
        
        # vectors for the solutions
        x_1 = np.zeros((len(t)))
        y_1 = np.zeros((len(t)))
        z_1 = np.zeros((len(t)))
        
        x_2 = np.zeros((len(t)))
        y_2 = np.zeros((len(t)))
        z_2 = np.zeros((len(t)))

        x_3 = np.zeros((len(t)))
        y_3 = np.zeros((len(t)))
        z_3 = np.zeros((len(t)))       

        x_4 = np.zeros((len(t)))
        y_4 = np.zeros((len(t)))
        z_4 = np.zeros((len(t)))
        
        x_5 = np.zeros((len(t)))
        y_5 = np.zeros((len(t)))
        z_5 = np.zeros((len(t)))

        x_6 = np.zeros((len(t)))
        y_6 = np.zeros((len(t)))
        z_6 = np.zeros((len(t)))
 
        # initial conditions
        """
        This initial conditions are random
        """
        #Define initial position vectors
        r1=[np.random.uniform(-2,2),np.random.uniform(-2,2),np.random.uniform(-2,2)]
        r2=[np.random.uniform(-2,2),np.random.uniform(-2,2),np.random.uniform(-2,2)]
        r3=[np.random.uniform(-2,2),np.random.uniform(-2,2),np.random.uniform(-2,2)]
        Pos_0=np.array([r1,r2,r3])
        #Convert pos vectors to arrays
        r1=np.array(r1,dtype="float64")
        r2=np.array(r2,dtype="float64")
        r3=np.array(r3,dtype="float64")
        #Define initial velocities
        v1=[np.random.uniform(-1,1),np.random.uniform(-1,1),np.random.uniform(-1,1)]
        v2=[np.random.uniform(-1,1),np.random.uniform(-1,1),np.random.uniform(-1,1)]
        v3=[np.random.uniform(-1,1),np.random.uniform(-1,1),np.random.uniform(-1,1)]

        #Convert velocity vectors to arrays
        v1=np.array(v1,dtype="float64")
        v2=np.array(v2,dtype="float64")
        v3=np.array(v3,dtype="float64")
        
        #Update COM formula
        r_com=(m1*r1+m2*r2+m3*r3)/(m1+m2+m3)
        #Update velocity of COM formula
        v_com=(m1*v1+m2*v2+m3*v3)/(m1+m2+m3)
        
        # transformation to CM
        v1 = v1-v_com
        v2 = v2-v_com
        v3 = v3-v_com
        Vel_0=np.array([v1-v_com,v2-v_com,v3-v_com])

        init_params=np.array([r1,r2,r3,v1,v2,v3]) #Initial parameters
        init_params=init_params.flatten() #Flatten to make 1D array

        three_body_sol=solve_ivp(ThreeBodyEquations, [0, 30], init_params, t_eval = t, method = 'Radau', max_step = 1e-2, args=(G,m1,m2,m3))

        r1_sol=three_body_sol.y[:3]
        r2_sol=three_body_sol.y[3:6]
        r3_sol=three_body_sol.y[6:9]
        
        x_1 = r1_sol[0]
        y_1 = r1_sol[1]
        z_1 = r1_sol[2]

        x_2 = r2_sol[0]
        y_2 = r2_sol[1]
        z_2 = r2_sol[2]
        
        x_3 = r3_sol[0]
        y_3 = r3_sol[1]
        z_3 = r3_sol[2]
        
       # Post Newtonian Solution

        u0 = np.zeros(18)


        u0[:3] = Pos_0[0,:]
        u0[3:6] = Vel_0[0,:]
        u0[6:9] = Pos_0[1,:]
        u0[9:12] = Vel_0[1,:]
        u0[12:15] = Pos_0[2,:]
        u0[15:] = Vel_0[2, :]
        
        three_body_sol=solve_ivp(dydx, [0, 30], u0, t_eval = t, method = 'Radau', max_step = 1e-2, args=(m1,m2,m3))

        r4_sol=three_body_sol.y[:3]
        r5_sol=three_body_sol.y[6:9]
        r6_sol=three_body_sol.y[12:15]
        
        x_4 = r4_sol[0]
        y_4 = r4_sol[1]
        z_4 = r4_sol[2]

        x_5 = r5_sol[0]
        y_5 = r5_sol[1]
        z_5 = r5_sol[2]
        
        x_6 = r6_sol[0]
        y_6 = r6_sol[1]
        z_6 = r6_sol[2]
        
        plot_solution(self, fig, x_1, y_1, z_1, t, x_2, y_2, z_2, x_3, y_3, z_3, 
         x_4, y_4, z_4, x_5, y_5, z_5, x_6, y_6, z_6)
           
        
if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    w = Widget()
    w.show()
    #main()
    sys.exit(app.exec_())
