import math
import csv
import os
from CoolProp.CoolProp import PropsSI as PSI
from scipy import constants

class Placas_aquecimento:

    def __init__(self, geometria, fluido1, fluido2, csv_path):
        self.geometria = geometria
        self.fluido1 = fluido1
        self.fluido2 = fluido2
        self.csv_path = csv_path
        self.parametros_placa()
        self.canais()
        self.parametros_portas()
        self.data = self.ler_csv()
        
    def canais(self):
        placas = self.geometria.get('placas', None)
        if placas is not None:
            self.canais_arrefecimento = int(placas / 2)
            self.canais_aquecimento = placas - self.canais_arrefecimento - 1
        else:
            raise ValueError ("O valor do número de placas não está corretamente definido na geometria do permutador.")
        
    def parametros_placa(self):
        d_g = self.geometria.get('d_g', None)
        plate_thickness = self.geometria.get('plate_thickness', None)
        plate_conductivity = self.geometria.get('plate_conductivity', None)
        l_v = self.geometria.get('l_v', None)
        l_p = self.geometria.get('l_p', None)
        w_p = self.geometria.get('w_p', None)
        l_h = self.geometria.get('l_h', None)
        plate_angle = self.geometria.get('plate_angle', None)
        plate_angle_radians = math.radians(plate_angle)
        wave_length = self.geometria.get('wave_length', None)
        area_catalgo = self.geometria.get('area_catalogo', None)

        if d_g is None or not isinstance(d_g, float):
            raise ValueError ("O espaçamento entre placas não está corretamente definido na geometria do permutador.")
        if plate_thickness is None or not isinstance(plate_thickness, float):
            raise ValueError ('A espessura da placa não está corretamente definida na geometria do permutador.')
        if plate_conductivity is None or not isinstance(plate_conductivity, float):
            raise ValueError ('A condutividade do material das placas não está corretamnete definida na geometria do permutador.')
        if l_v is None or not isinstance(l_v, float):
            raise ValueError ('A altura do permutador não está corretamente definida na sua geometria.')
        if area_catalgo is None or not isinstance(area_catalgo, float):
            raise ValueError ('A área do permutador não está corretamente definida na sua geometria.')
        if l_p is None or not isinstance(l_p, float):
            raise ValueError ('A distância vertical entre diâmetros das portas de fluído não está corretamente definida na geometria do permutador.')
        if w_p is None or not isinstance(w_p, float):
            raise ValueError ('A largura do permutador não está corretamente definida na sua geometria.')
        if l_h is None or not isinstance(l_h, float):
            raise ValueError ('A distância horizontal entre diâmetros das portas de fluído não está corretamente definida na geometria do permutador.')
        if wave_length is None or not isinstance(wave_length, float):
            raise ValueError ('O comprimento de onda da corrugação da placa não está corretamente definido na geometria do permutador.')
        
        self.compressed_pitch = d_g + plate_thickness
        self.amplitude = d_g / 2
        omega = (math.pi * d_g) / wave_length
        enlargement_factor = (1 / 6) * (1 + math.sqrt(1 + (omega ** 2)) + 4 * math.sqrt(1 + ((omega ** 2) / 2)))
        self.d_h = 2 * d_g / enlargement_factor
        self.d_e = 2 * d_g
        self.corrugation_aspect_ratio = 2 * d_g / wave_length
        self.area_cross = d_g * w_p
        
        return {
            "d_g": d_g,
            "plate_thickness": plate_thickness,
            "l_v": l_v,
            "plate_conductivity": plate_conductivity,
            "l_p": l_p,
            "w_p": w_p,
            "l_h": l_h,
            "plate_angle": plate_angle,
            "plate_angle_radians": plate_angle_radians,
            "wave_length": wave_length,
            "area_catalogo": area_catalgo,
        }        
        
    def parametros_portas(self):
        port_D_wf = self.geometria.get('port_D_wf', None)
        port_D_cooling = self.geometria.get('port_D_cooling', None)
        if port_D_wf is None or not isinstance(port_D_wf, float):
            raise ValueError ('O diâmetro da porta de entrada do fluído refrigerante não está corretamente definido.')
        if port_D_cooling is None or not isinstance(port_D_cooling, float):
            raise ValueError ('O diâmetro da porta de saída do fluído de arrefecimento não está corretamente definido.')

        return {
            "port_D_wf": port_D_wf,
            "port_D_cooling": port_D_cooling,
        }
    
    def ler_csv(self):
        data = {}

        with open(self.csv_path, newline='', encoding='utf-8') as csvfile:
            reader = csv.reader(csvfile, delimiter=';')
            headers = next(reader)
            filtered_headers = [header for header in headers if "unit" not in header]
            expected_vars = ['m', 'v', 'p', 'h', 'T', 'Td_bp', 'vol', 'x', 's']
            var_indices = {var: None for var in expected_vars}
            for index, header in enumerate(filtered_headers):
                if header in var_indices:
                    var_indices[header] = index
            if any(index is None for index in var_indices.values()):
                raise ValueError("Falta de informação no CSV.")
            for row in reader:
                row_data = {}
                for var, index in var_indices.items():
                    value = row[index]  # Get the raw value first
                    try:
                        row_data[var] = float(value)  # Try to convert to float
                    except ValueError:
                        # print(f"Warning: Could not convert '{value}' to float for variable '{var}'. Keeping it as string.")
                        row_data[var] = value  # Keep it as string
                label = row[len(filtered_headers) - 1]
                data[label] = row_data

        return data

    def valor_ligacao(self, label, var_name):
        if label in self.data:
            return self.data[label].get(var_name)
        else:
            raise ValueError(f"Label '{label}' não foi encontrado no CSV.")
        
    def calculo_areas_perdas(self, c1, c2, c10, c11):
        c1_data = self.data.get(c1.label)
        c10_data = self.data.get(c10.label)
        c2_data = self.data.get(c2.label)
        c11_data = self.data.get(c11.label)

        q1 = c1_data['m'] * (c1_data['h'] * 1000 - PSI("H", "Q", 1, "P", c1_data['p'], self.fluido1))
        T_cooling_mid1 = c10_data['T'] + q1 / (c10_data['m'] * PSI("C", "T", T_cooling_mid1 + 273.15, "P", c10_data['p'], self.fluido2))

        c_fluido1_1 = c1_data['m'] * PSI("C", "T", 273.15 + ((PSI("T", "Q", 1, "P", c1_data['p'], self.fluido1) - 273.15) + c1_data['T']) / 2, "P", c1_data['p'], self.fluido1)
        c_fluido2_1 = c10_data['m'] * PSI("C", "T", 273.15 + (T_cooling_mid1 + c10_data['T']) / 2, "P", c10_data['p'], self.fluido2)
        c_min_1 = min(c_fluido1_1, c_fluido2_1)
        c_max_1 = max(c_fluido1_1, c_fluido2_1)
        Q_1_max = c_min_1 * (c1_data['T'] - c10_data['T'])
        efficiency_1 = q1 / Q_1_max
        Cr_1 = c_min_1 / c_max_1
        ntu_1 = (1 / (Cr_1 - 1)) * math.log((efficiency_1 - 1) / (efficiency_1 * Cr_1 - 1))

        if c2_data['x'] == -1:
            q2 = c1_data['m'] * (PSI("H", "Q", 1, "P", c1_data['p'], self.fluido1) - PSI("H", "Q", 0, "P", c1_data['p'], self.fluido1))
            T_cooling_mid2 = T_cooling_mid1 + q2 / (c10_data['m'] * PSI("C", "T", T_cooling_mid1 + 273.15, "P", c10_data['p'], self.fluido2))
            c_water_2 = c10_data['m'] * PSI("C", "T", 273.15 + (T_cooling_mid1 + T_cooling_mid2) / 2, "P", c10_data['p'], self.fluido2)
            Q_2_max = c_water_2 * ((PSI("T", "Q", 1, "P", c1_data['p'], self.fluido1) - 273.15) - T_cooling_mid1)
            efficiency_2 = q2 / Q_2_max
            ntu_2 = -math.log(1 - efficiency_2)

            q3 = abs(c1_data['Q']) - q2 - q1
            c_wf_3 = c1_data['m'] * PSI("C", "T", 273.15 + ((PSI("T", "Q", 0, "P", c1_data['p'], self.fluido1) - 273.15) + c2_data['T']) / 2, "P", c1_data['p'], self.fluido1)  # W
            c_water_3 = c10_data['m'] * PSI("C", "T", 273.15 + (T_cooling_mid2 + c11_data['T']) / 2, "P", c10_data['p'], self.fluido2)  # W
            c_min_3 = min(c_wf_3, c_water_3)
            c_max_3 = max(c_wf_3, c_water_3)
            Cr_3 = c_min_3 / c_max_3
            Q_3_max = c_min_3 * ((PSI("T", "Q", 0, "P", c1_data['p'], self.fluido1) - T_cooling_mid2 - 273.15))
            efficiency_3 = q3 / Q_3_max
            ntu_3 = (1 / (Cr_3 - 1)) * math.log((efficiency_3 - 1) / (efficiency_3 * Cr_3 - 1))
        else:
            q2 = c1_data['m'] * (PSI("H", "Q", 1, "P", c1_data['p'], self.fluido1) - PSI("H", "Q", 0, "P", c1_data['p'], self.fluido1)) * (1 - c2_data['x'])
            T_cooling_mid2 = c11_data['T']
            c_water_2 = c10_data['m'] * PSI("C", "T", 273.15 + (T_cooling_mid1 + T_cooling_mid2) / 2, "P", c10_data['p'], self.fluido2)  # W
            Q_2_max = c_water_2 * ((PSI("T", "Q", 1, "P", c1_data['p'], self.fluido1) - 273.15) - T_cooling_mid1)
            efficiency_2 = q2 / Q_2_max
            ntu_2 = -math.log(1 - efficiency_2)
            q3 = 0
            ntu_3 = 0

        A_flow_wf = self.area_cross * self.canais_arrefecimento
        A_flow_water = self.area_cross * self.canais_aquecimento

        G_fluido1 = c1_data['m'] / A_flow_wf
        G_fluido2 = c10_data['m'] / A_flow_water

        Re_canal_fluido1_1 = G_fluido1 * self.d_h / PSI("V", "T", 273.15 + (c1_data['T'] + PSI("T", "P", c1_data['p'], "Q", 1, self.fluido1) - 273.15) / 2, "P", c1_data['p'], self.fluido1)

        if c2_data['x'] == -1:
            Re_canal_fluido1_3 = G_fluido1 * self.d_h / PSI("V", "T", 273.15 + (PSI("T", "P", c1_data['p'], "Q", 0, self.fluido1) + c2_data['T'] - 273.15) / 2, "P", c1_data['p'], self.fluido1)
            Re_canal_fluido2_3 = G_fluido2 * self.d_h / PSI("V", "T", 273.15 + (T_cooling_mid2 + c11_data['T']) / 2, "P", c10_data['p'], self.fluido2)
        else:
            Re_canal_fluido1_3 = 0
            Re_canal_fluido2_3 = 0

        Re_canal_fluido2_1 = G_fluido2 * self.d_h / PSI("V", "T", 273.15 + (c10_data['T'] + T_cooling_mid1) / 2, "P", c10_data['p'], self.fluido2)
        Re_canal_fluido2_2 = G_fluido2 * self.d_h / PSI("V", "T", 273.15 + (T_cooling_mid1 + T_cooling_mid2) / 2, "P", c10_data['p'], self.fluido2)

        R_f_h = 0.0001
        R_f_c = 0.0001

        Pr_fluido2_1 = (PSI("V", "T", 273.15 + (c10_data['T'] + c11_data['T']) / 2, "P", c10_data['p'], self.fluido2) * PSI("C", "T", 273.15 + (c10_data['T'] + c11_data['T']) / 2, "P", c10_data['p'], self.fluido2) / PSI("L", "T", 273.15 + (c10_data['T'] + c11_data['T']) / 2, "P", c10_data['p'], self.fluido2))
        Pr_fluido1_1 = PSI("V", "T", 273.15 + (c1_data['T'] + PSI("T", "P", c1_data['p'], "Q", 1, self.fluido1) - 273.15) / 2, "P", c1_data['p'], self.fluido1) * PSI("C", "T", 273.15 + (c1_data['T']+ PSI("T", "P", c1_data['p'], "Q", 1, self.fluido1) - 273.15) / 2, "P", c1_data['p'], self.fluido1) / PSI("L", "T", 273.15 + (c1_data['T'] + PSI("T", "P", c1_data['p'], "Q", 1, self.fluido1) - 273.15) / 2, "P", c1_data['p'], self.fluido1)
        h_fluido2_1 = 0.724 * (PSI("L", "T", 273.15 + (c10_data['T'] + c11_data['T']) / 2, "P", c10_data['p'], self.fluido2) / self.d_h) * ((self.plate_angle / 30) ** (0.646)) * (Re_canal_fluido2_1 ** 0.583) * (Pr_fluido2_1 ** (1/3))
        h_fluido1_1 = 0.724 * (PSI("L", "T", 273.15 + (c1_data['T'] + PSI("T", "P", c10_data['p'], "Q", 1, self.fluido1) - 273.15) / 2, "P", c1_data['p'], self.fluido1) / self.d_h) * ((self.plate_angle / 30) ** (0.646)) * (Re_canal_fluido1_1 ** 0.583) * (Pr_fluido1_1 ** (1/3))
        u_1 = 1 / ((1 / h_fluido2_1) + (1 / h_fluido1_1) + (self.area_catalogo * R_f_h) + (self.area_catalogo * R_f_c))
        a_1 = ntu_1 * c_min_1 / u_1

        if c2_data['x'] == -1:
            titulo_medio = 0.5
        else:
            titulo_medio = (1 + c2_data['x']) / 2

        racio_densidade_fluido1 = PSI("D", "Q", 0, "P", c1_data['p'], self.fluido1) / PSI("D", "Q", 1, "P", c1_data['p'], self.fluido1)    
        Pr_fluido2_2 = PSI("V", "T", 273.15 + (T_cooling_mid1 + T_cooling_mid2) / 2, "P", c10_data['p'], self.fluido2) * PSI("C", "T", 273.15 + (T_cooling_mid1 + T_cooling_mid2) / 2, "P", c10_data['p'], self.fluido2) / PSI("L", "T", 273.15 + (T_cooling_mid1 + T_cooling_mid2) / 2, "P", c10_data['p'], self.fluido2)
        h_fluido2_2 = 0.724 * (PSI("L", "T", 273.15 + (T_cooling_mid1 + T_cooling_mid2) / 2, "P", c10_data['p'], self.fluido2) / self.d_h) * ((self.plate_angle / 30) ** (0.646)) * (Re_canal_fluido2_2 ** 0.583) * (Pr_fluido2_2 ** (1/3))

        Bd = (PSI("D", "Q", 0, "P", c1_data['p'], self.fluido1) - PSI("D", "Q", 1, "P", c1_data['p'], self.fluido1)) * constants.g * (self.d_h ** 2) / PSI("I", "Q", 0, "P", c1_data['p'], self.fluido1)
        Pr_l = PSI("C", "Q", 0, "P", c1_data['p'], self.fluido1) * PSI("V", "Q", 0, "P", c1_data['p'], self.fluido1) / PSI("L", "Q", 0, "P", c1_data['p'], self.fluido1)
        densidade_fluido1_media = 1 / ((titulo_medio / PSI("D", "Q", 1, "P", c1_data['p'], self.fluido1)) + (1 - titulo_medio) / PSI("D", "Q", 0, "P", c1_data['p'], self.fluido1))
        G_fluido1_eq = G_fluido1 * (1 - titulo_medio + titulo_medio * (racio_densidade_fluido1 ** 0.5))
        Re_fluid1_eq = G_fluido1_eq * self.d_h / PSI("V", "Q", 0, "P", c1_data['p'], self.fluido1)
        Nu_fluido1_bi_fase = 4.3375 * (Re_fluid1_eq ** 0.5383) ** (Pr_l ** (1/3)) * (Bd ** (-0.3872))
        h_fluido1_bi_fase = Nu_fluido1_bi_fase * (PSI("L", "Q", 0, "P", c1_data['p'], self.fluido1)) / self.d_h

        Pr_fluido2_2 = PSI("V", "T", 273.15 + (T_cooling_mid1 + T_cooling_mid2) / 2, "P", c10_data['p'], self.fluido2) * PSI("C", "T", 273.15 + (T_cooling_mid1 + T_cooling_mid2) / 2, "P", c10_data['p'], self.fluido2) / PSI("L", "T", 273.15 + (T_cooling_mid1 + T_cooling_mid2) / 2, "P", c10_data['p'], self.fluido2)
        h_fluido2_2 = 0.724 * (PSI("L", "T", 273.15 + (T_cooling_mid1 + T_cooling_mid2) / 2, "P", c10_data['p'], self.fluido2) / self.d_h) * ((self.plate_angle / 30) ** (0.646)) * (Re_canal_fluido2_2 ** 0.583) * (Pr_fluido2_2 ** (1/3))
        u_2 = 1 / ((1 / h_fluido2_2) + (1 / h_fluido1_bi_fase) + (self.area_catalogo * R_f_h) + (self.area_catalogo * R_f_c))
        a_2 = ntu_2 * c_water_2 / u_2

        if c2_data['x'] == -1:
            Pr_fluido1_3 = PSI("V", "T", 273.15 + (PSI("T", "P", c1_data['p'], "Q", 0, self.fluido1) - 273.15 + c2_data['T']) / 2, "P", c1_data['p'], self.fluido1) * PSI("C", "T", 273.15 + (PSI("T", "P", c1_data['p'], "Q", 0, self.fluido1) - 273.15 + c2_data['T']) / 2, "P", c1_data['p'], self.fluido1) / PSI("L", "T", 273.15 + (PSI("T", "P", c1_data['p'], "Q", 0, self.fluido1) - 273.15 + c2_data['T']) / 2, "P", c1_data['p'], self.fluido1)
            h_fluido1_3 = 0.724 * (PSI("L", "T", 273.15 + (PSI("T", "P", c1_data['p'], "Q", 0, self.fluido1) - 273.15 + c2_data['T']) / 2, "P", c1_data['p'], self.fluido1) / self.d_h) * ((self.plate_angle / 30) ** (0.646)) * (Re_canal_fluido1_3 ** 0.583) * (Pr_fluido1_3 ** (1/3))
            Pr_fluido2_3 = PSI("V", "T", 273.15 + (T_cooling_mid2 + c11_data['T']) / 2, "P", c10_data['p'], self.fluido2) * PSI("C", "T", 273.15 + (T_cooling_mid2 + c11_data['T']) / 2, "P", c10_data['p'], self.fluido2) / PSI("L", "T", 273.15 + (T_cooling_mid2 + c11_data['T']) / 2, "P", c10_data['p'], self.fluido2)
            h_fluido2_3 = 0.724 * (PSI("L", "T", 273.15 + (T_cooling_mid2 + c11_data['T']) / 2, "P", c10_data['p'], self.fluido2) / self.d_h) * ((self.plate_angle / 30) ** (0.646)) * (Re_canal_fluido2_3 ** 0.583) * (Pr_fluido2_3 ** (1/3))
            u_3 = 1 / ((1 / h_fluido2_3) + (1 / h_fluido1_3) + (self.area_catalogo * R_f_h) + (self.area_catalogo * R_f_c))
            a_3 = ntu_3 * c_min_3 / u_3
        else:
            a_3 = 0

        area_needed = a_1 + a_2 + a_3

        We = (G_fluido1 ** 2) * self.d_h / (densidade_fluido1_media * PSI("I", "P", c1_data['p'], "Q", 0, self.fluido1))

        if Re_canal_fluido2_1 < 2000:
            f_0 = 16 / Re_canal_fluido2_1
            f_1 = (149 / Re_canal_fluido2_1) + 0.9625
        else:
            f_0 = (1.56 * math.log(Re_canal_fluido2_1) - 3) ** (-2)
            f_1 = 9.75 / (Re_canal_fluido2_1 ** 0.289)

        f_fluido2_1 = ((math.cos(math.radians(90 - self.plate_angle))) / (0.045 * math.tan(math.radians(90 - self.plate_angle)) + 0.09 * math.sin(math.radians(90 - self.plate_angle)) + f_0 / math.cos(math.radians(90 - self.plate_angle))) + ((1 - math.cos(math.radians(90 - self.plate_angle))) / (math.sqrt(3.8 * f_1)))) ** (-2)

        if Re_canal_fluido2_2 < 2000:
            f_0 = 16 / Re_canal_fluido2_2
            f_1 = (149 / Re_canal_fluido2_2) + 0.9625
        else:
            f_0 = (1.56 * math.log(Re_canal_fluido2_2) - 3) ** (-2)
            f_1 = 9.75 / (Re_canal_fluido2_2 ** 0.289)

        f_fluido2_2 = ((math.cos(math.radians(90 - self.plate_angle))) / (0.045 * math.tan(math.radians(90 - self.plate_angle)) + 0.09 * math.sin(math.radians(90 - self.plate_angle)) + f_0 / math.cos(math.radians(90 - self.plate_angle))) + ((1 - math.cos(math.radians(90 - self.plate_angle))) / (math.sqrt(3.8 * f_1)))) ** (-2)
        
        if c2_data['x'] == -1:
            if Re_canal_fluido2_3 < 2000:
                f_0 = 16 / Re_canal_fluido2_3
                f_1 = (149 / Re_canal_fluido2_3) + 0.9625
            else:
                f_0 = (1.56 * math.log(Re_canal_fluido2_3) - 3) ** (-2)
                f_1 = 9.75 / (Re_canal_fluido2_3 ** 0.289)
            
            f_fluido2_3 = ((math.cos(math.radians(90 - self.plate_angle))) / (0.045 * math.tan(math.radians(90 - self.plate_angle)) + 0.09 * math.sin(math.radians(90 - self.plate_angle)) + f_0 / math.cos(math.radians(90 - self.plate_angle))) + ((1 - math.cos(math.radians(90 - self.plate_angle))) / (math.sqrt(3.8 * f_1)))) ** (-2)
        else:
            f_0 = 0
            f_1 = 0
        f_fluido2_3 = 0

        if Re_canal_fluido1_1 < 2000:
            f_0 = 16 / Re_canal_fluido1_1
            f_1 = (149 / Re_canal_fluido1_1) + 0.9625
        else:
            f_0 = (1.56 * math.log(Re_canal_fluido1_1) - 3) ** (-2)
            f_1 = 9.75 / (Re_canal_fluido1_1 ** 0.289)
        
        f_fluido1_1 = ((math.cos(math.radians(90 - self.plate_angle))) / (0.045 * math.tan(math.radians(90 - self.plate_angle)) + 0.09 * math.sin(math.radians(90 - self.plate_angle)) + f_0 / math.cos(math.radians(90 - self.plate_angle))) + ((1 - math.cos(math.radians(90 - self.plate_angle))) / (math.sqrt(3.8 * f_1)))) ** (-2)
        f_fluido1_bi_fase = 0.0146 * (Re_fluid1_eq ** 0.9814) * (We ** (-1.0064))    

        if c2_data['x'] == -1:
            if Re_canal_fluido1_3 < 2000:
                f_0 = 16 / Re_canal_fluido1_3
                f_1 = (149 / Re_canal_fluido1_3) + 0.9625
            else:
                f_0 = (1.56 * math.log(Re_canal_fluido1_3) - 3) ** (-2)
                f_1 = 9.75 / (Re_canal_fluido1_3 ** 0.289)

            f_fluido1_3 = ((math.cos(math.radians(90 - self.plate_angle))) / (0.045 * math.tan(math.radians(90 - self.plate_angle)) + 0.09 * math.sin(math.radians(90 - self.plate_angle)) + f_0 / math.cos(math.radians(90 - self.plate_angle))) + ((1 - math.cos(math.radians(90 - self.plate_angle))) / (math.sqrt(3.8 * f_1)))) ** (-2)   
        else:
            f_0 = 0
            f_1 = 0
        f_fluido1_3 = 0

        delta_p_f_fluido2_1 = (4 * f_fluido2_1 * self.l_p * (G_fluido2 ** 2)) / (2 * PSI("D", "P", c10_data['p'], "T", 273.15 + T_cooling_mid1, self.fluido2) * self.d_h)
        delta_p_f_fluido2_2 = (4 * f_fluido2_2 * self.l_p * (G_fluido2 ** 2)) / (2 * PSI("D", "P", c10_data['p'], "T", 273.15 + T_cooling_mid2, self.fluido2) * self.d_h)
        if c2_data['x'] == -1:
            delta_p_f_fluido2_3 = (4 * f_fluido2_3 * self.l_p * (G_fluido2 ** 2)) / (2 * PSI("D", "P", c10_data['p'], "T", 273.15 + (c11_data['T'] + T_cooling_mid2) / 2, self.fluido2) * self.d_h)
        else:
            delta_p_f_fluido2_3 = 0

        if c2_data['x'] == -1:
            delta_p_f_fluido1_3 = (4 * f_fluido1_3 * self.l_p * (G_fluido1 ** 2)) / (2 * PSI("D", "P", c1_data['p'], "T", 273.15 + (c2.T.val + PSI("T", "P", c1_data['p'], "Q", 0, self.fluido1) - 273.15) / 2, self.fluido1) * self.d_h)
        else:
            delta_p_f_fluido1_3 = 0

        delta_p_f_fluido1_2 = (4 * f_fluido1_bi_fase * self.l_p * (G_fluido1 ** 2)) / (2 * ((PSI("D", "P", c1_data['p'], "Q", 1, self.fluido1) + PSI("D", "P", c1_data['p'], "Q", 0, self.fluido1)) / 2) * self.d_h)
        delta_p_f_fluido1_1 = (4 * f_fluido1_1 * self.l_p * (G_fluido1 ** 2)) / (2 * PSI("D", "P", c1_data['p'], "T", 273.15 + (c1.T.val + PSI("T", "Q", 1, "P", c1_data['p'], self.fluido1) - 273.15) / 2, self.fluido1) * self.d_h)

        delta_p_f_fluido1 = delta_p_f_fluido1_1 + delta_p_f_fluido1_2 + delta_p_f_fluido1_3
        delta_p_f_fluido2 = delta_p_f_fluido2_1 + delta_p_f_fluido2_2 + delta_p_f_fluido2_3

        G_porta_fluido1 = 4 * c1_data['m'] / (math.pi * (self.port_D_wf ** 2))
        G_porta_fluido2 = 4 * c10_data['m'] / (math.pi * (self.port_D_cooling ** 2))

        delta_p_p_fluido2 = 1.4 * (G_porta_fluido2 ** 2) * 2 / (2 * PSI("D", "P", c10_data['p'], "T", 273.15 + (c10_data['T'] + c11_data['T']) / 2, self.fluido2))
        delta_p_p_fluido1 = 1.4 * (G_porta_fluido1 ** 2) * 2 / ((PSI("D", "P", c1_data['p'], "Q", 1, self.fluido1) + PSI("D", "P", c1_data['p'], "Q", 0, self.fluido1)))

        delta_p_fluido2 = delta_p_f_fluido2 + delta_p_p_fluido2
        delta_p_fluido1 = delta_p_f_fluido1 + delta_p_p_fluido1
        
        self.delta_p_fluido1 = delta_p_fluido1
        self.delta_p_fluido2 = delta_p_fluido2

        return {
            area_needed,
            delta_p_fluido2,
            delta_p_fluido1,
        }
        
    
        
    
        



        
        
        
        
