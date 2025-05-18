import numpy as np
import matplotlib.pyplot as plt
def read_initial_conditions(filename):
    # ux_i uy_i vx_i vy_i pour chaque nœud
    return np.loadtxt(filename)

def analytical_solution(u0, v0, omega, t):
    """Retourne (u(t), v(t)) pour une vibration harmonique."""
    u = u0 * np.cos(omega * t) + v0 / omega * np.sin(omega * t)
    v = -u0 * omega * np.sin(omega * t) + v0 * np.cos(omega * t)
    return u, v

def read_numerical_solution(filename):
    # t_i ux_i uy_i vx_i vy_i pour un seul noeud
    return np.loadtxt(filename)

def compute_truncation_error(numerical, initial_conditions, omega):
    times = numerical[:, 0]
    h= numerical[1, 0] - numerical[0, 0]  # Pas de temps
    u0x, u0y, v0x, v0y = initial_conditions
    errors = []

    for i  in range(len(times)-1):
        t, ux_num, uy_num, vx_num, vy_num = numerical[i]
        t1, ux_num1, uy_num1, vx_num1, vy_num1 = numerical[i+1]

        ux_exact, vx_exact = analytical_solution(u0x, v0x, omega, t)
        uy_exact, vy_exact = analytical_solution(u0y, v0y, omega, t)
        ux_exact1,_= analytical_solution(u0x, v0x, omega, t + h)
     
     

        err = {
            "time": t,
            "err_ux": abs(ux_num - ux_exact),
            "err_uy": abs(uy_num - uy_exact),
            "err_vx": abs(vx_num - vx_exact),
            "err_vy": abs(vy_num - vy_exact),
            # "err_t_x": abs((ux_exact1 - ux_exact)/h -vx_exact), 
            "err": np.sqrt((ux_exact-ux_num)**2 + (uy_exact-uy_num)**2),
            
        }
        errors.append(err)

    return errors

def main():
    omega = 1665.0  # fréquence naturelle (ex. 10 rad/s)
    
    initial_conditions = read_initial_conditions("initial_fork_1.0.txt")[0]  # On prend un seul nœud
    numerical = read_numerical_solution("time_fork.txt")

    errors = compute_truncation_error(numerical, initial_conditions, omega)

    times = [e['time'] for e in errors]
    err_ux = [e['err_ux'] for e in errors]
    err_uy = [e['err_uy'] for e in errors]
    err_vx = [e['err_vx'] for e in errors]
    err_vy = [e['err_vy'] for e in errors]
    # err_t = [e['err_t_x'] for e in errors]
    errr= [e['err'] for e in errors]

    for e in errors:
        print(f"t={e['time']:.5f}", f"err_ux={e['err_ux']:.3e}, err_uy={e['err_uy']:.3e}", 
              f"err_vx={e['err_vx']:.3e}, err_vy={e['err_vy']:.3e}")
    plt.plot(times, err_ux, label="Erreur ux")
    plt.plot(times, err_uy, label="Erreur uy")
    # plt.plot(times, np.ones(len(times))*np.mean(err_ux), label="Erreur vx")
    # print("Erreur ux moyenne:", np.mean(err_ux))
    plt.plot(times, err_vx, label="Erreur vx")
    # plt.plot(times, errr, label="Erreur totale")
    plt.plot(times, err_vy, label="Erreur vy")
    # plt.plot(times,err_t, label="Erreur troncature", linestyle='-')
    # plt.plot(times,np.mean(errr)*np.ones(len(times)), label="Erreur moyenne vx", linestyle='--')
    plt.xlabel("Temps [s]")
   
  
    plt.ylabel("Erreur absolue")
    plt.title("Erreurs en fonction du temps")
    plt.legend()
    plt.grid()
    plt.savefig("erreurs.pdf")  # Décommentez pour sauvegarder la figure
    plt.show()
if __name__ == "__main__":
    main()