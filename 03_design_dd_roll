import pybert as pb
import pygimli as pg
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata


def dataformatsyscal(data, name='suveyfile.txt'):
    """function to convert a data file into a format readable by electepro"""
    import math
    # f = open('surveyfile.txt', 'w')
    f = open(name, 'w')

    f.write('#  X   Y   Z \r\n')

    for i in range(data.sensorCount()):
        x = '%.2f' % (data.sensorPositions()[i][0])
        y = '%.2f' % (data.sensorPositions()[i][1])
        z = '%.2f' % (data.sensorPositions()[i][2])

        f.write(str(i + 1) + '  ' + str(x) + '  ' + str(y) + '   ' + str(z) + '\r\n')

    f.write('#  A   B   M   N \r\n')

    for i in range(len(data('a'))):
        f.write(str(i + 1) + '  ' + str(int(data('a')[i]) + 1) + '  ' + str(int(data('b')[i]) + 1) + '  ' + str(
            int(data('m')[i]) + 1) + '  ' + str(int(data('n')[i]) + 1) + '\r\n')

    f.close()

def pseudo_section_with_depth(data, other_array=False, array=None, array_type='dd', interp_method='linear', scatter=True, cmin=None, cmax=None):
    '''function to plot pseudo sections at depth
    depth of investigation was obtained from Loke practical guide
    array_type can be :
    - wa = wenner alpha
    - dd
    interp_method can be:
    - nearest
    - cubic
    - linear'''
    nconfig = len(data('a'))
    pos = np.array(data.sensorPositions())
    xe_list = []
    ze_list = []
    rhoa_list = []
    try:
        if other_array:
            rhoa_from_data = array
        else:
            rhoa_from_data = data.get('rhoa').array()
        scatter_only = False
        if all(res_v == 0 for res_v in rhoa_from_data):
            scatter_only = True
    except:
        print('No apparent resistivity values found')
        rhoa_from_data = data.get('a').array() * 0
        scatter_only = True

    for i in range(nconfig):
        a = int(data('a')[i])
        b = int(data('b')[i])
        m = int(data('m')[i])
        n = int(data('n')[i])

        x_c1 = pos[a, 0]
        y_c1 = pos[a, 1]
        x_c2 = pos[b, 0]
        y_c2 = pos[b, 1]
        x_p1 = pos[m, 0]
        y_p1 = pos[m, 1]
        x_p2 = pos[n, 0]
        y_p2 = pos[n, 1]

        if array_type == 'dd':
            a_dist = x_c2 - x_c1
            n_dist = (x_p1 - x_c2) / a_dist
            n_vals = np.array([1, 2, 3, 4, 5, 6, 7, 8])
            ze_over_a_vals = np.array([0.416, 0.697, 0.962, 1.220, 1.476, 1.730, 1.983, 2.236])
            ze_over_a = np.interp(n_dist, n_vals, ze_over_a_vals)
            ze = -ze_over_a * a_dist
            xe = (x_c1 + x_c2 + x_p1 + x_p2) / 4
            rhoa = rhoa_from_data[i]
        elif array_type == 'wa':
            a_dist = x_p1 - x_c1
            ze_over_a = 0.519
            ze = -ze_over_a * a_dist
            xe = (x_c1 + x_c2 + x_p1 + x_p2) / 4
            rhoa = rhoa_from_data[i]

        else:
            print('Array cannot be recognized')

        xe_list.append(xe)
        ze_list.append(ze)
        rhoa_list.append(rhoa)

    # Contour plot
    xe_array = np.array(xe_list).reshape(len(xe_list), 1)
    ze_array = np.array(ze_list).reshape(len(ze_list), 1)
    xze = np.hstack((xe_array, ze_array))
    rhoa_array = np.array(rhoa_list)
    npoints = 1000
    # npoints_x = int(np.abs((np.max(pos[:, 0]) - np.min(pos[:, 1]))) / (pos[1, 0] - pos[0, 0]))
    # npoints_z = int(np.abs((np.max(ze_array)) - np.min(ze_array)) / (pos[1, 0] - pos[0, 0]))

    if scatter_only:
        plt.figure()
        plt.scatter(xe_array, ze_array, c='k', alpha=1, marker='.')
        plt.title('Scatter plot')

    else:
        xi = np.linspace(np.min(xe_array), np.max(xe_array), npoints).reshape(npoints, 1)
        zi = np.linspace(np.min(ze_array), np.max(ze_array), npoints).reshape(npoints, 1)
        X, Z = np.meshgrid(xi, zi)

        grid_rhoa = griddata(xze, rhoa_array, (X, Z), method=interp_method)

        v = np.linspace(np.min(rhoa_array), np.max(rhoa_array), 100, endpoint=True)
        plt.figure()
        CS = plt.contourf(X, Z, grid_rhoa, v, vmin=cmin, vmax=cmax, cmap='jet')
        plt.colorbar(CS)
        plt.scatter(xe_array, ze_array, c='k', alpha=0.1, marker='.')
        plt.title('pseudo section ' + array_type)

        if scatter:
            plt.figure()
            plt.plot(rhoa_array, 'x')
            plt.title('Apparent resistivity')

def geomfactor(data, index, printfactor=True):
    a = int(data('a')[index])
    b = int(data('b')[index])
    m = int(data('m')[index])
    n = int(data('n')[index])

    xa = data.sensorPositions()[a][0]
    ya = data.sensorPositions()[a][1]
    za = data.sensorPositions()[a][2]

    xb = data.sensorPositions()[b][0]
    yb = data.sensorPositions()[b][1]
    zb = data.sensorPositions()[b][2]

    xm = data.sensorPositions()[m][0]
    ym = data.sensorPositions()[m][1]
    zm = data.sensorPositions()[m][2]

    xn = data.sensorPositions()[n][0]
    yn = data.sensorPositions()[n][1]
    zn = data.sensorPositions()[n][2]

    ram = ((xa - xm) ** 2 + (ya - ym) ** 2 + (za - zm) ** 2) ** 0.5
    ran = ((xa - xn) ** 2 + (ya - yn) ** 2 + (za - zn) ** 2) ** 0.5
    rbm = ((xb - xm) ** 2 + (yb - ym) ** 2 + (zb - zm) ** 2) ** 0.5
    rbn = ((xb - xn) ** 2 + (yb - yn) ** 2 + (zb - zn) ** 2) ** 0.5

    rapm = ((xa - xm) ** 2 + (-ya - ym) ** 2 + (za - zm) ** 2) ** 0.5
    rapn = ((xa - xn) ** 2 + (-ya - yn) ** 2 + (za - zn) ** 2) ** 0.5
    rbpm = ((xb - xm) ** 2 + (-yb - ym) ** 2 + (zb - zm) ** 2) ** 0.5
    rbpn = ((xb - xn) ** 2 + (-yb - yn) ** 2 + (zb - zn) ** 2) ** 0.5

    H = 1 / ram - 1 / ran - 1 / rbm + 1 / rbn + 1 / rapm - 1 / rapn - 1 / rbpm + 1 / rbpn
    if H == 0.0:
        K = 0.0
    else:
        K = 4 * np.pi / H

    if printfactor:
        print('the geometric factor is:', K)
    return K

def dd_generation(scheme, max_geom_factor):
    '''function to generate dipole-dipole scheme with varying dipole spacing
    geometric factor is approximated with linear topography'''

    ne = scheme.sensorCount()
    scheme.resize(1)
    a_list = []
    b_list = []
    m_list = []
    n_list = []

    db = 1
    while db < ne:
        dbm = db
        for i in range(ne):
            dbm = db
            for j in range(ne):
                a = i
                b = i + db
                m = i + db + dbm
                n = i + db + dbm + db

                scheme('a')[0] = a
                scheme('b')[0] = b
                scheme('m')[0] = m
                scheme('n')[0] = n
                try:
                    geo_factor = geomfactor(scheme, 0, printfactor=False)
                except:
                    geo_factor = max_geom_factor + 1


                if (np.abs(geo_factor) < max_geom_factor) and (n < ne):
                    a_list.append(a)
                    b_list.append(b)
                    m_list.append(m)
                    n_list.append(n)

                else:
                    break

                # dbm = dbm + 1 # the separation is to short
                dbm = dbm + db

        db = db + 1

    scheme.resize(len(a_list))
    for i in range(len(a_list)):
        scheme('a')[i] = a_list[i]
        scheme('b')[i] = b_list[i]
        scheme('m')[i] = m_list[i]
        scheme('n')[i] = n_list[i]
        scheme('valid')[i] = 1

# Reference file
x_ini = 0.0
end_line = 143
ne = 144
scheme_ref = pb.createData(elecs=pg.utils.grange(start=x_ini, end=end_line, n=ne), schemeName=None)
dd_generation(scheme_ref, max_geom_factor=5000)
dataformatsyscal(scheme_ref, name='dd_reference.txt')
nconfig_ref = scheme_ref.size()

# Roll-along file
x_ini_roll = 36
end_line_roll = 179
ne = 144
scheme_roll = pb.createData(elecs=pg.utils.grange(start=x_ini_roll, end=end_line_roll, n=ne), schemeName=None)
dd_generation(scheme_roll, max_geom_factor=5000)
nconfig_roll = scheme_roll.size()
pseudo_section_with_depth(scheme_ref, array_type='dd')

# Remove overlap configurations
for iroll in range(nconfig_roll):
    print(iroll, ' configurations analyzed out of: ', nconfig_roll)
    xa_roll = scheme_roll.sensorPositions()[int(scheme_roll.get('a')[iroll])][0]
    xb_roll = scheme_roll.sensorPositions()[int(scheme_roll.get('b')[iroll])][0]
    xm_roll = scheme_roll.sensorPositions()[int(scheme_roll.get('m')[iroll])][0]
    xn_roll = scheme_roll.sensorPositions()[int(scheme_roll.get('n')[iroll])][0]
    for iref in range(nconfig_ref):
        xa_ref = scheme_ref.sensorPositions()[int(scheme_ref.get('a')[iref])][0]
        xb_ref = scheme_ref.sensorPositions()[int(scheme_ref.get('b')[iref])][0]
        xm_ref = scheme_ref.sensorPositions()[int(scheme_ref.get('m')[iref])][0]
        xn_ref = scheme_ref.sensorPositions()[int(scheme_ref.get('n')[iref])][0]

        if xa_roll == xa_ref and xb_roll == xb_ref and xm_roll == xm_ref and xn_roll == xn_ref:
            scheme_roll('valid')[iroll] = 0
            break
            
scheme_roll.removeInvalid()
pseudo_section_with_depth(scheme_roll, array_type='dd')
dataformatsyscal(scheme_roll, name='dd_roll.txt')

plt.show()



