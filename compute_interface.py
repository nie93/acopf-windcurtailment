import math


def calculatePg(f):
    if not f:
        return None
    
    device_id = list(f)

    vg1_m = f[device_id[0]]['VAYPM_mag'] / 7.5
    vg2_m = f[device_id[0]]['VBYPM_mag'] / 7.5
    vg3_m = f[device_id[0]]['VAZPM_mag'] / 37.5
    vg4_m = f[device_id[0]]['VBZPM_mag'] / 37.5
    ig1_m = f[device_id[0]]['IAWPM_mag'] * 10. / 3
    ig2_m = f[device_id[0]]['IBWPM_mag'] * 8. / 15
    ig3_m = f[device_id[0]]['ICWPM_mag'] * 10.
    ig4_m = f[device_id[0]]['IAXPM_mag'] * 20. / 3

    vg1_a = f[device_id[0]]['VAYPM_ang']
    vg2_a = f[device_id[0]]['VBYPM_ang']
    vg3_a = f[device_id[0]]['VAZPM_ang']
    vg4_a = f[device_id[0]]['VBZPM_ang']
    ig1_a = f[device_id[0]]['IAWPM_ang']
    ig2_a = f[device_id[0]]['IBWPM_ang']
    ig3_a = f[device_id[0]]['ICWPM_ang']
    ig4_a = f[device_id[0]]['IAXPM_ang']
    

    pg1 = 3 * vg1_m * ig1_m * math.cos(vg1_a - ig1_a)
    pg2 = 3 * vg2_m * ig2_m * math.cos(vg2_a - ig2_a)
    pg3 = 3 * vg3_m * ig3_m * math.cos(vg3_a - ig3_a)
    pg4 = 3 * vg4_m * ig4_m * math.cos(vg4_a - ig4_a)

    return [pg1, pg2, pg3, pg4]
    