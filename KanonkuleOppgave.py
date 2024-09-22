# Importerer libraries
import numpy as np
import matplotlib.pyplot as plt

# Utregnings funksjon, v0 er farten vi setter inn, valg er om vi vil ha luftmotstand med
def simulate(v0, valg):
# Konstanter som angår personen/simulasjonen
    m = 75  # Masse i kg
    g = 9.81  # Gravitasjonell akselerasjon i m/s^2
    D = 1 # Diameter i m
    b = 1.6*10**(-4)*D  # Lineær dragkoeffisient
    c = 0.25*D**2  # Kvadratiske dragkoeffisienten
    tmax = 100  # Maks tid i sekunder
    dt = 0.01  # Tids steg hver runthrough
    t = 0
# Starts variabler som angår situasjonen vi er i ved start
    theta = np.radians(45)  # Avskyter vinkel i radianer
    vx = v0 * np.cos(theta) # Starts fart i x retning
    vy = v0 * np.sin(theta) # Starts fart i y retning
    x = 0 # Starts posisjon langs x aksen
    y = 5 * np.sqrt(2) # Starts høyde
    px = vx * m # Starts bevegelses mengde i x retning
    py = vy * m # Starts bevegels esmengde i y retning
# Tomme lister for plotting
    x_vals = []
    y_vals = []
# Selve simulasjons loopen
# Stoper når vi går mer enn tmax sekunder eller når vi faller under bakken
    while t < tmax and y > 0:
        t += dt # Øker tidssteg
        # Oppdaterer fart
        vx = px / m
        vy = py / m
        # Oppdaterer krefter
        Fx = (-b * vx - c * abs(vx) * vx) * valg
        Fy = -m * g + (-b * vy - c * abs(vy) * vy) * valg
        # Oppdaterer bevegelses mengde
        px += Fx * dt
        py += Fy * dt
        # Oppdaterer posisjonene
        x += vx * dt
        y += vy * dt
        # Legger posisjonene i listene for plotting
        x_vals.append(x)
        y_vals.append(y)
    return x_vals, y_vals # Etter loopen er ferdig returner listene


# Newtonian funksjon for å finne riktig starts verdi
def newtonian(valg):
    goal = 50  # Målets avstand i meter
    step = 20  # Starts steg lengde, justeres hvis man skyter over/under
    v0 = 20  # Vilkårlig gjettet starts hastighet
    flip = True # Når målet passeres flippes denne og minker steg lengden
    hit = False # Loopen går til målet blir truffet
# Selve loopen
    while not hit:
        x_vals, y_vals = simulate(v0, valg) # Får posisjonslistene fra simulasjonen
        if valg==1: # Hvis vi har med luftmotstand legges farten i en liste for plotting
            v_list1.append(v0)
        elif valg==0: # Dersom vi ikke har det med legges det i en annen liste for plotting
            v_list2.append(v0)
        if round(x_vals[-1], 1) < goal: # Dersom vi treffer under målet
            if not flip: # Hvis flip er false gjør steglengde mindre
                step *= 0.5
                flip = True # Og sett til true
            v0 += step # Øker v0 med step
            print(f"Under given: {v0}") # Printer ut at vi er når under målet
        elif round(x_vals[-1], 1) > goal: # Treffer vi over målet
            if flip: # Hvis flip er True gjør steglengden mindre
                step *= 0.5
                flip = False # Sett flip til False
            v0 -= step # Minker v0 med step
            print(f"Over given: {v0}") # Printer ut av vi traff over målet
# Dersom forskjellen mellom der vi lander og målet er mindre enn 4cm, eller steglengden blir for liten
        if abs(x_vals[-1] - goal) <= 0.04 or step <= 0.001:
            print(f"HUZZZAA! Correct velocity: {v0} m/s") # Print ut verdien vi fikk
            print(f"Final x position: {x_vals[-1]} meters") # Og hvor vi traff
            hit = True # Bryter while løkken
    return v0 # Returnerer verdien vi fant

# Plotter
fig, ax = plt.subplots(2, 2, figsize=(10, 10))
ax = ax.flatten()

# Liste med hvordan farten utviklet seg over tid
v_list1=[]
v_list2=[]

# Løser for fart vi trenger med luftmotstand
print("-Solving for with drag-")
v_drag = newtonian(1) # valg=1 så vi har med luftmotstand
# Simulerer hva som skjer hvis vi bruker farten fi fant hvis vi bruker/ignorerer luftmotstand
xw_drag, yw_drag = simulate(v_drag, 1) # Med
xn_drag, yn_drag = simulate(v_drag, 0) # Uten

# Plotter resultatene på samme graf for å vise forskjellen
ax[0].plot(xw_drag, yw_drag, "b.", label="Med luftmotstand")
ax[0].plot(xn_drag, yn_drag, "r", label="Uten luftmotstand")
ax[0].set_title(f"Med luftmotstand trenger vi v0={round(v_drag,1)}m/s")
ax[0].set_xlabel("x (m)")
ax[0].set_ylabel("y (m)")
ax[0].grid()
ax[0].legend()


# Løser for farten vi trenger uten luftmotstand
print("\n-Solving for without drag-")
v_nodrag = newtonian(0) # Valg=0 så vi har ikke med luftmotstand
# Simulerer hva som skjer hvis vi bruker farten vi fant med/uten luftmotstand
x2w_drag, y2w_drag = simulate(v_nodrag, 1) # Med
x2n_drag, y2n_drag = simulate(v_nodrag, 0) # Uten

# Plotter resultatne på samme graf for å vise forskjellen
ax[1].plot(x2w_drag, y2w_drag, "b", label="Med luftmotstand")
ax[1].plot(x2n_drag, y2n_drag, "r.", label="Uten luftmotstand")
ax[1].set_title(f"Uten luftmotstand trenger vi v0={round(v_nodrag,1)}m/s")
ax[1].set_xlabel("x (m)")
ax[1].set_ylabel("y (m)")
ax[1].grid()
ax[1].legend()

# Plotter hvordan farten utviklet seg over tid når vi løste for med/uten luftmotstand
ax[2].plot(v_list1,"green")
ax[2].set_title("Finner fart v0 med luftmotstand")
ax[2].set_xlabel("Forsøk")
ax[2].set_ylabel("v0")
ax[2].grid()
ax[2].legend()

ax[3].plot(v_list2,"green")
ax[3].set_title("Finner fart v0 uten luftmotstand")
ax[3].set_xlabel("Forsøk")
ax[3].set_ylabel("v0")
ax[3].grid()
ax[3].legend()

plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.show()

# Printer ut resultatene
print(f"-When solving for with drag-\nVelocity needed when drag accounted for:  {round(v_drag,1)}m/s\nDistance reached when drag accounted for: {round(xw_drag[-1],5)}m\nDistance reached when drag ignored:       {round(xn_drag[-1],5)}m\n")

print(f"-When solving for without drag-\nVelocity needed when drag is ignored:     {round(v_nodrag,1)}m/s\nDistance reached when drag accounted for: {round(x2w_drag[-1],5)}m\nDistance reached when drag ignored:       {round(x2n_drag[-1],5)}m")


