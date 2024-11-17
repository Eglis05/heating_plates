#Foram comentadas linhas do source code em src/tespy/components/heat_exchangers/base.py para remoção das mensagens caso não dê para calcular o valor da eficiência térmica (linhas 1073 e 1086, versão 0.7.6.post1)

from tespy.networks import Network
from tespy.tools.characteristics import CharLine
from CoolProp.CoolProp import PropsSI as PSI
import matplotlib.pyplot as plt
from fluprodia import FluidPropertyDiagram
import math
import os
import numpy as np
from scipy import constants
from scipy.special import iv, kv
from tespy.connections import Connection
from tespy.components import (CycleCloser, Compressor, Valve, Pump, HeatExchanger, Source, Sink)
from permutadores import Placas_aquecimento as placas

heatPump = Network(T_unit='C', p_unit='Pa', h_unit='kJ / kg', iterinfo=False)

closer = CycleCloser('Cycle closer') 
Cooling1 = HeatExchanger('Permutador de placas Heating / Radiador Cooling')
Cooling2 = HeatExchanger('Radiador Heating / Permutador de placas Cooling')
safety = Valve('Válvula segurança')
needle1 = Valve('Válvula agulha 1')
needle2 = Valve('Válvula agulha 2')
valve1 = Valve('Expansion valve 1') 
valve2 = Valve('Expansion valve 2') 
valve_three_wayA = Valve('Válvula 3 vias A')
valve_three_wayB = Valve('Válvula 3 vias B')
circuit_filter = Valve('Filtro')
compressor = Compressor('Compressor') 
pump = Pump('Pump') 

source1 = Source('ambIn1')
source2 = Source('ambIn2')
sink1 = Sink('ambOut1')
sink2 = Sink('ambOut2')

c0 = Connection(circuit_filter, 'out1', compressor, 'in1', label= 'Estado inicial')
c1 = Connection(compressor, 'out1', safety, 'in1', label= 'Pós compressão')
c02 = Connection(safety, 'out1', valve_three_wayA, 'in1', label= 'Válvula segurança')
c03 = Connection(valve_three_wayA, 'out1', needle1, 'in1', label= 'Saída válvula 3 vias')
c04 = Connection(needle1, 'out1', Cooling1, 'in1', label= 'Válvula agulha 1')
c2 = Connection(Cooling1, "out1", closer, "in1", label= "Pós arrefecimento")
c3 = Connection(closer, "out1", needle2, "in1", label= "Cycle Closer")
c05 = Connection(needle2, 'out1', valve1, 'in1', label= 'Válvula agulha 2')
c31 = Connection(valve1, "out1", valve2, "in1", label= "Pós expansão inicial")
c4 = Connection(valve2, "out1", Cooling2, "in2", label= "Pós expansão final")
c00 = Connection(Cooling2, 'out2', valve_three_wayB, 'in1', label= 'Saída permutador ar')
c01 = Connection(valve_three_wayB, 'out1', circuit_filter, 'in1', label= 'Saída válvula 3 vias / Entrada no filtro')

c10 = Connection(source1, "out1", Cooling1, "in2", label= "Água")
c11 = Connection(Cooling1, "out2", sink1, "in1", label= "Água aquecida")

c20 = Connection(source2, "out1", Cooling2, "in1", label= "Ambiente antes")
c21 = Connection(Cooling2, "out1", sink2, "in1", label= "Ambiente depois")

heatPump.add_conns(c0, c1, c02, c03, c04, c2, c3, c05, c31, c4, c00, c01, c10, c11, c20, c21)

wf = "R32"
massRateWF = 0.037
T_wf_in = -25
pressureIn = PSI("P", "T", 273.15 - 26, "Q", 1, wf)
pressureOut = PSI("P", "T", 273.15 + 15, "Q", 0, wf)
pressureRatio = pressureOut / pressureIn
pressureMix = 4e5
h0 = PSI("H", "P", pressureIn, "T", 273.15 + T_wf_in, wf) / 1000
fluidCheck = PSI("T", "P", pressureIn, "Q", 1, wf) - 273.15
T_amb = -15
T_cooling_in1 = 5
T_cooling_in2 = T_amb
coolingHeated1 = 12
mix = 'INCOMP::GKN[0.6]'
coolingFluid1 = mix
coolingFluid2 = "Air"
coolingVolumeRate1 = 2.8 / 3600 #m3/s
coolingVolumeRate2 = 5200 / 3600 #m3/s
coolingDensity1 = PSI("D", "T", 273.15 + T_cooling_in1, "P", pressureMix, coolingFluid1) #kg/m3
massRateCooling1 = coolingVolumeRate1 * coolingDensity1 #kg/s
massRateCooling2 = PSI("D", "T", 273.15 + T_amb, "P", 1.01325e5, "Air") * coolingVolumeRate2

kv_safety = 5.5
kv_needle1 = 0.85 * 0.865
kv_needle2 = 0.85 * 0.865
delta_p_filter = 0.21e5
delta_p_3way = 0.1e5

t_pack = T_amb
t_pack_final = 20
m_pack = 335
cp_pack = 500
a_pack = 1.073 #m
b_pack = 0.578 #m
c_pack = 0.586 #m
volume_pack = a_pack * b_pack * c_pack #m3
transfer_area_pack = 4 * a_pack * c_pack + a_pack * b_pack #m2

if T_wf_in <= fluidCheck:
    print("Existência de líquido / mistura bifásica no compressor. Aumente a temperatura de entrada ou reduza a pressão inicial.")
else:
    print("Modo de aquecimento")
    condenserPressureLossRatioHot = 1
    evaporatorPressureLossRatioHot = 1

    efficiencyCompressor = 0.85

    efficiencyPump = 0.8

    Cooling1.set_attr(pr1 = 1, pr2 = 1)
    compressor.set_attr(eta_s = efficiencyCompressor)
    Cooling2.set_attr(pr1 = 1, pr2 = 1)

    c0.set_attr(fluid={wf: 1}, T = T_wf_in, p = pressureIn, m = massRateWF)
    c2.set_attr(x = 0, T = 20)
    c31.set_attr(x = 0.15)
    safety.set_attr(pr = 1)
    valve_three_wayA.set_attr(pr = 1)
    valve_three_wayB.set_attr(pr = 1)    
    needle1.set_attr(pr = 1)
    needle2.set_attr(pr = 1)
    circuit_filter.set_attr(pr = 1)

    c10.set_attr(fluid={coolingFluid1: 1}, p = pressureMix, T = T_cooling_in1)
    c11.set_attr(T = coolingHeated1)

    c20.set_attr(fluid={coolingFluid2: 1}, p = 1.01325e5, T = T_cooling_in2)
    c21.set_attr(m = massRateCooling2)
    
    h =  PSI("H", "T", 273.15 + T_cooling_in1, "P", pressureMix, coolingFluid1) / 1e3
    [c.set_attr(h0=h0) for c in [c10, c11]]

    heatPump.solve(mode='design')

    density_safety = PSI("D", "P", c1.p.val, "T", 273.15 + c1.T.val, wf)
    safety.set_attr(pr = 1 - (kv_safety * (3600 * massRateWF / density_safety)) / c1.p.val)

    heatPump.solve(mode='design')
    valve_three_wayA.set_attr(pr = (c02.p.val - delta_p_3way) / c02.p.val)

    heatPump.solve(mode='design') 

    density_needle1 = PSI("D", "P", c03.p.val, "T", 273.15 + c03.T.val, wf)
    needle1.set_attr(pr = 1 - (kv_needle1 * (3600 * massRateWF / density_needle1)) / c03.p.val)

    heatPump.solve(mode='design') 
    
    density_needle2 = PSI("D", "P", c2.p.val, "Q", c2.x.val, wf)
    needle2.set_attr(pr = 1 - (kv_needle2 * (3600 * massRateWF / density_needle2)) / c2.p.val)

    heatPump.solve(mode='design') 

    valve_three_wayB.set_attr(pr = (c00.p.val - delta_p_3way) / c00.p.val)

    heatPump.solve(mode='design') 

    circuit_filter.set_attr(pr = (c01.p.val - delta_p_filter) / c01.p.val)

    heatPump.solve(mode='design')

    c0.set_attr(p = None)
    c4.set_attr(p = pressureIn)

    heatPump.solve(mode='design')
    heatPump.save('hp_design_heating')
    data = os.path.join('hp_design_heating', 'connections.csv')

    cop = abs(Cooling1.Q.val) / compressor.P.val

    #Dimensionamento do permutador de placas paralelas B3-027

    plates = 30
    area_catalogo = (plates - 2) * 0.026

    geometria = {
        "port_D_wf": 0.0127,          
        "port_D_cooling": 0.0127,   
        "plate_thickness": 0.0005,    
        "plate_conductivity": 14.9,  
        "plate_angle": 60,
        "wave_length": 0.008,        
        "l_v": 0.311,           
        "l_p": 0.25,           
        "w_p": 0.111, 
        "l_h": 0.05,                   
        "d_g": 0.00236,                
        "placas":  plates,
        "area_catalogo": area_catalogo,                   
    }

    permutador = placas(geometria, wf, coolingFluid1, data)
    permutador.calculo_areas_perdas(c0, c1, c10, c11)

    delta_p_c_ratio = 1 - (permutador.delta_p_fluido2 / c10.p.val)
    delta_p_h_ratio = 1 - (permutador.delta_p_fluido1 / c2.p.val)

    while abs((Cooling1.pr2.val - delta_p_c_ratio) / Cooling1.pr2.val) > 0.00001:
        Cooling1.pr2.val = (Cooling1.pr2.val + delta_p_c_ratio) / 2
        heatPump.solve(mode='design')

    while abs((Cooling1.pr1.val - delta_p_h_ratio) / Cooling1.pr1.val) > 0.00001:
        Cooling1.pr1.val = (Cooling1.pr1.val + delta_p_h_ratio) / 2
        heatPump.solve(mode='design')

    heatPump.save('hp_design_heating')
    heatPump.print_results()

    result_dict = {}    
    result_dict.update({Cooling2.label: Cooling2.get_plotting_data()[2]})
    result_dict.update({compressor.label: compressor.get_plotting_data()[1]})
    result_dict.update({Cooling1.label: Cooling1.get_plotting_data()[1]})
    result_dict.update({valve1.label: valve1.get_plotting_data()[1]})
    result_dict.update({valve2.label: valve2.get_plotting_data()[1]})
    result_dict.update({valve_three_wayA: valve_three_wayA.get_plotting_data()[1]})
    result_dict.update({valve_three_wayB: valve_three_wayB.get_plotting_data()[1]})
    result_dict.update({needle1.label: needle1.get_plotting_data()[1]})
    result_dict.update({needle2.label: needle2.get_plotting_data()[1]})
    result_dict.update({circuit_filter.label: circuit_filter.get_plotting_data()[1]})
    result_dict.update({safety.label: safety.get_plotting_data()[1]})

    diagram = FluidPropertyDiagram(wf)
    diagram.set_unit_system(T='°C', p='Pa', h='kJ/kg')

    for key, data in result_dict.items():
        result_dict[key]['datapoints'] = diagram.calc_individual_isoline(**data)

    diagram.calc_isolines()

    fig, ax = plt.subplots(1, figsize=(16, 10))
    diagram.draw_isolines(fig, ax, 'logph', x_min= c3.h.val - 350, x_max= c1.h.val + 100, y_min= 1e4  , y_max= c1.p.val * 10)

    for key in result_dict.keys():
        datapoints = result_dict[key]['datapoints']
        ax.plot(datapoints['h'], datapoints['p'], color='#ff0000')
        ax.scatter(datapoints['h'][0], datapoints['p'][0], color='#ff0000')
    fig.savefig(str(wf) + '_logph_heating.svg')

    print('O COP da bomba de calor com ' + str(wf)+ ' é ' + str(round(cop, 3)))
    print('O rácio de compressão no aquecimento é ' + str(round(pressureRatio, 2)))
    print("No aquecimento, a velocidade do fluído de trabalho à saída do compressor é " + str(round(4 * massRateWF / (math.pi * PSI("D", "T", 273.15 + c1.T.val, "P", pressureOut, wf) * (0.0127 ** 2)), 2)) + " m/s")
    print("Área necessária (permutador de placas): " + str(round(permutador.area_needed, 3)) + " m2")
    print("Área catalogo (permutador de placas): " + str(round(area_catalogo,3)) + " m2")
    print("Perda de carga no lado da água: " + str(round(permutador.delta_p_fluido2 / 1000, 4)) + " kPa")
    print("Rácio de pressão no lado da água: " + str(round(delta_p_c_ratio, 3)))
    print("Perda de carga no lado do fluído de trabalho: " + str(round(permutador.delta_p_fluido1 / 1000, 3)) + " kPa")
    print("Rácio de pressão no lado da fluído de trabalho após permutador de placas: " + str(round(delta_p_h_ratio, 3)))

    water_rate_heating = str(round(c10.m.val / PSI("D", "T", 273.15 + (c10.T.val + c11.T.val) / 2, "P", pressureMix, mix) * 3600, 2))
    water_rate_heating_v = str(round((c10.m.val / PSI("D", "T", 273.15 + (c10.T.val + c11.T.val) / 2, "P", pressureMix, mix) * 3600) / (0.25 * (geometria.get.port_D_cooling ** 2) * math.pi * 3600), 2))

    # print("\n")

    # print("Modo arrefecimento")

    # heatPump_cooling= Network(T_unit='C', p_unit='Pa', h_unit='kJ / kg', iterinfo=False)

    # #Dimensionamento do permutador de placas paralelas B3-018

    # port_D_wf = 0.0127 #m
    # port_D_cooling = 0.0127 #m
    # plate_thickness = 0.0005 #m / d_p
    # plate_conductivity = 14.9 # W/m.K
    # chevron_angle = 60 #
    # chevron_angle_radians = math.radians(chevron_angle)
    # l_v = 0.231 #m
    # l_p = 0.182 #m
    # # l_p = l_v - port_D #m
    # w_p = 0.09 #m
    # l_h = 0.043 #m
    # d_g = 0.00226 # Espaço médio entre placas / channel gap / b #m
    # # corrugation_pitch = 0.0075 #
    # plates = 28
    # channels_cooling = int(plates / 2)
    # channels_refrigerant = plates - channels_cooling - 1
    # area_catalogo = 0.018 * (plates - 2) #m2
    # compressed_pitch = d_g + plate_thickness #m
    # amplitude = d_g / 2 #m
    # wave_length = 0.008 #m experimentação
    # omega = (math.pi * d_g) / wave_length #m
    # enlargement_factor = (1 / 6) * (1 + math.sqrt(1 + (omega ** 2)) + 4 * math.sqrt(1 + ((omega ** 2) / 2)))
    # d_h = 2 * d_g / enlargement_factor #m
    # d_e = 2 * d_g #m
    # corrugation_aspect_ratio = 2 * d_g / wave_length #m
    # area_cross = d_g * w_p #m2

    # T_cooling_in1 = 26
    # T_wf_in = 12
    # T_amb = 60
    # pressureIn = PSI("P", "T", 273.15 + 10, "Q", 1, wf)
    # pressureOut = PSI("P", "T", 273.15 + 63, "Q", 1, wf)
    # pressureRatio_cooling = pressureOut / pressureIn
    # coolingHeated1 = 14
    # coolingVolumeRate1 = 3.6 / 3600 #m3/s    
    # coolingVolumeRate2 = 5200 / 3600 #m3/s
    # massRateCooling1 = PSI("D", "T", 273.15 + T_cooling_in1, "P", pressureMix, mix) * coolingVolumeRate1    
    # coolingDensity2 = PSI("D", "T", 273.15 + T_amb, "P", 1.01325e5, "Air") #kg/m3 
    # massRateCooling2 = coolingVolumeRate2 * coolingDensity2 #kg/s

    # c0 = Connection(circuit_filter, 'out1', compressor, 'in1', label= 'Estado inicial')
    # c1 = Connection(compressor, 'out1', safety, 'in1', label= 'Pós compressão')
    # c02 = Connection(safety, 'out1', valve_three_wayA, 'in1', label= 'Válvula segurança')
    # c03 = Connection(valve_three_wayA, 'out1', Cooling1, 'in1', label= 'Saída válvula 3 vias')
    # c2 = Connection(Cooling1, "out1", closer, "in1", label= "Pós arrefecimento")
    # c3 = Connection(closer, "out1", valve1, "in1", label= "Cycle Closer")
    # c31 = Connection(valve1, "out1", valve2, "in1", label= "Pós expansão inicial")
    # c4 = Connection(valve2, "out1", needle1, "in1", label= "Pós expansão final")
    # c04 = Connection(needle1, 'out1', Cooling2, 'in2', label= 'Válvula agulha 1')
    # c05 = Connection(Cooling2, 'out2', needle2, 'in1', label= 'Pós permutador de placas')
    # c06 = Connection(needle2, 'out1', valve_three_wayB, 'in1', label= 'Válvula de agulha 2')
    # c07 = Connection(valve_three_wayB, 'out1', circuit_filter, 'in1', label= 'Saída válvula 3 vias / Entrada no filtro')

    # c10 = Connection(source1, "out1", Cooling1, "in2", label= "Ambiente antes")
    # c11 = Connection(Cooling1, "out2", sink1, "in1", label= "Ambiente depois")

    # c20 = Connection(source2, "out1", Cooling2, "in1", label= "Água antes")
    # c21 = Connection(Cooling2, "out1", sink2, "in1", label= "Água depois")

    # heatPump_cooling.add_conns(c0, c1, c02, c03, c2, c3, c31, c4, c04, c05, c06, c07, c10, c11, c20, c21)

    # Cooling1.set_attr(pr1 = 1, pr2 = 1)
    # Cooling2.set_attr(pr1 = 1, pr2 = 1)
    # # compressor.set_attr(pr = pressureRatio_cooling)
    # circuit_filter.set_attr(pr = 1)
    # safety.set_attr(pr = 1)
    # needle1.set_attr(pr = 1)
    # needle2.set_attr(pr = 1)
    # valve_three_wayA.set_attr(pr = 1)
    # valve_three_wayB.set_attr(pr = 1)

    # c0.set_attr(fluid={wf: 1}, T = T_wf_in, p = pressureIn, m = massRateWF)
    # c2.set_attr(x = 0, T = 63)
    # c31.set_attr(x = 0.2)

    # c10.set_attr(fluid={coolingFluid2: 1, coolingFluid1: None}, p = 1.01325e5, T = T_amb, m = massRateCooling2)
    # # c11.set_attr(T = 30)
    # c20.set_attr(fluid={coolingFluid1: 1, coolingFluid2: None}, p = pressureMix, T = T_cooling_in1)
    # c21.set_attr(T = coolingHeated1)

    # h =  PSI("H", "T", 273.15 + T_cooling_in1, "P", pressureMix, coolingFluid1) / 1e3
    # [c.set_attr(h0=h0) for c in [c20, c21]]

    # heatPump_cooling.solve(mode= 'design')

    # density_safety = PSI("D", "P", c1.p.val, "T", 273.15 + c1.T.val, wf)
    # safety.set_attr(pr = 1 - (kv_safety * (3600 * massRateWF / density_safety)) / pressureOut)

    # heatPump_cooling.solve(mode= 'design')

    # valve_three_wayA.set_attr(pr = (c02.p.val - delta_p_3way) / c02.p.val)    

    # heatPump_cooling.solve(mode= 'design')

    # density_needle1 = PSI("D", "P", c4.p.val, "Q", c4.x.val, wf)
    # needle1.set_attr(pr = 1 - (2 * kv_needle1 * (3600 * massRateWF / density_needle1)) / c4.p.val)

    # heatPump_cooling.solve(mode= 'design')

    # density_needle2 = PSI("D", "P", c05.p.val, "T", 273.15 + c05.T.val, wf)
    # needle2.set_attr(pr = 1 - (2 * kv_needle2 * (3600 * massRateWF / density_needle1)) / c05.p.val)

    # heatPump_cooling.solve(mode='design')

    # valve_three_wayB.set_attr(pr = (c06.p.val - delta_p_3way) / c06.p.val)

    # heatPump_cooling.solve(mode='design')

    # circuit_filter.set_attr(pr = (c07.p.val - delta_p_filter) / c07.p.val) 

    # heatPump_cooling.solve(mode='design') 

    # c0.set_attr(p = None)
    # c4.set_attr(p = pressureIn)

    # eer = abs(Cooling2.Q.val) / compressor.P.val

    # if c4.x.val == -1:
    #     average_quality_in_out = 0.5
    #     q1 = c0.m.val * (PSI("H", "Q", 0, "P", c1.p.val, wf) - c4.h.val * 1000) #W
    #     T_cooling_mid1 = c20.T.val - q1 / (c20.m.val * PSI("C", "T", c20.T.val + 273.15, "P", c20.p.val, coolingFluid1))
    #     c_wf_1 = c0.m.val * PSI("C", "T", (c1.T.val + PSI("T", "Q", 1, "P", pressureIn, wf) - 273.15) / 2 + 273.15, "P", pressureIn, wf) #W/K
    #     c_water_1 = c20.m.val * PSI("C", "T", 273.15 + (T_cooling_mid1 + c20.T.val) / 2, "P", c20.p.val, coolingFluid1) #W/K
    #     c_min_1 = min(c_wf_1, c_water_1)
    #     c_max_1 = max(c_wf_1, c_water_1)
    #     Q_1_max = c_min_1 * (c20.T.val - c4.T.val)
    #     efficiency_1 = q1 / Q_1_max
    #     Cr_1 = c_min_1 / c_max_1
    #     ntu_1 = (1 / (Cr_1 - 1)) * math.log((efficiency_1 - 1) / (efficiency_1 * Cr_1 - 1))

    #     q2 = c0.m.val * (PSI("H", "Q", 1, "P", pressureIn, wf) - PSI("H", "Q", 0, "P", pressureIn, wf))
    #     T_cooling_mid2 = T_cooling_mid1 - q2 / (c20.m.val * PSI("C", "T", T_cooling_mid1 + 273.15, "P", c20.p.val, coolingFluid1))
    #     c_water_2 = c20.m.val * PSI("C", "T", 273.15 + (T_cooling_mid1 + T_cooling_mid2) / 2, "P", c20.p.val, coolingFluid1) #W
    #     Q_2_max = c_water_2 * (T_cooling_mid1 - (PSI("T", "Q", 1, "P", pressureIn, wf) - 273.15))
    #     efficiency_2 = q2 / Q_2_max
    #     ntu_2 = - math.log(1 - efficiency_2)

    # else:
    #     q1 = 0
    #     ntu_1 = 0
    #     c_min_1 = 0
    #     average_quality_in_out = (1 + c4.x.val) / 2
    #     q2 = c0.m.val * (PSI("H", "Q", 1, "P", pressureIn, wf) - c4.h.val * 1000)
    #     T_cooling_mid1 = c20.T.val
    #     T_cooling_mid2 = T_cooling_mid1 - q2 / (c20.m.val * PSI("C", "T", T_cooling_mid1 + 273.15, "P", c20.p.val, coolingFluid1))
    #     c_water_2 = c20.m.val * PSI("C", "T", 273.15 + (T_cooling_mid1 + T_cooling_mid2) / 2, "P", c20.p.val, coolingFluid1) #W
    #     Q_2_max = c_water_2 * (T_cooling_mid1 - c4.T.val)
    #     efficiency_2 = q2 / Q_2_max
    #     ntu_2 = - math.log(1 - efficiency_2)

    # q3 = c0.m.val * (c0.h.val * 1000 - PSI("H", "Q", 1, "P", pressureIn, wf)) #W
    # c_wf_3 = c0.m.val * PSI("C", "T", 273.15 + ((PSI("T", "Q", 0, "P", pressureIn, wf) - 273.15) + c0.T.val) / 2, "P", pressureIn, wf) #W
    # c_water_3 = c20.m.val * PSI("C", "T", 273.15 + (c21.T.val + T_cooling_mid2) / 2, "P", c20.p.val, coolingFluid2) #W
    # c_min_3 = min(c_wf_3, c_water_3)
    # c_max_3 = max(c_wf_3, c_water_3)    
    # Q_3_max = c_min_3 * (T_cooling_mid2 - (PSI("T", "Q", 1, "P", pressureIn, wf) - 273.15))
    # Cr_3 = c_min_3 / c_max_3    
    # efficiency_3 = q3 / Q_3_max
    # ntu_3 = (1 / (Cr_3 - 1)) * math.log((efficiency_3 - 1) / (efficiency_3 * Cr_3 - 1))

    # G_port_h = 4 * c0.m.val / (math.pi * (port_D_wf ** 2))
    # G_port_c = 4 * c20.m.val / (math.pi * (port_D_cooling ** 2))

    # A_flow_wf = area_cross * channels_refrigerant
    # A_flow_water = area_cross * channels_cooling

    # G_channel_wf = c0.m.val / A_flow_wf
    # G_channel_water = c20.m.val / A_flow_water
    # Re_channel_wf3 = G_channel_wf * d_h / PSI("V", "T", 273.15 + (PSI("T", "P", pressureIn, "Q", 1, wf) - 273.15 + c0.T.val) / 2, "P", pressureIn, wf)
    # Re_channel_water = G_channel_water * d_h / PSI("V", "T", 273.15 + (c20.T.val + c21.T.val) / 2, "P", c20.p.val, coolingFluid1)

    # if c4.x.val == -1:
    #     Pr_c1 = PSI("V", "T", 273.15 + (c20.T.val + c21.T.val) / 2, "P", c20.p.val, coolingFluid1) * PSI("C", "T", 273.15 + (c20.T.val + c21.T.val) / 2, "P", c20.p.val, coolingFluid1) / PSI("L", "T", 273.15 + (c20.T.val + c21.T.val) / 2, "P", c20.p.val, coolingFluid1)
    #     Pr_h1 = PSI("V", "T", 273.15 + (c4.T.val + PSI("T", "P", pressureIn, "Q", 0, wf) - 273.15 ) / 2, "P", pressureIn, wf) * PSI("C", "T", 273.15 + (c4.T.val + PSI("T", "P", pressureIn, "Q", 0, wf) - 273.15) / 2, "P", pressureIn, wf) / PSI("L", "T", 273.15 + (c4.T.val + PSI("T", "P", pressureIn, "Q", 0, wf) - 273.15) / 2, "P", pressureIn, wf)
    #     h_c1 = 0.724 * (PSI("L", "T", 273.15 + (c20.T.val + c21.T.val) / 2, "P", c20.p.val, coolingFluid1) / d_h) * ((chevron_angle / 30) ** (0.646)) * (Re_channel_water ** 0.583) * (Pr_c1 ** (1/3))
    #     h_1 = 0.724 * (PSI("L", "T", 273.15 + (c4.T.val + PSI("T", "P", pressureIn, "Q", 0, wf)) / 2, "P", pressureIn, wf) / d_h) * ((chevron_angle / 30) ** (0.646)) * (Re_channel_wf1 ** 0.583) * (Pr_h1 ** (1/3))
    #     u_1 = 1 / ((1 / h_c1) + (1 / h_1) + (area_catalogo * R_f_h) + (area_catalogo * R_f_c))
    #     a_1 = ntu_1 * c_min_1 / u_1
    # else:
    #     a_1 = 0

    # Pr_c2 = PSI("V", "T", 273.15 + (T_cooling_mid1 + T_cooling_mid2) / 2, "P", c20.p.val, coolingFluid1) * PSI("C", "T", 273.15 + (T_cooling_mid1 + T_cooling_mid2) / 2, "P", c20.p.val, coolingFluid1) / PSI("L", "T", 273.15 + (T_cooling_mid1 + T_cooling_mid2) / 2, "P", c20.p.val, coolingFluid1)
    # h_c2 = 0.724 * (PSI("L", "T", 273.15 + (T_cooling_mid1 + T_cooling_mid2) / 2, "P", c20.p.val, coolingFluid1) / d_h) * ((chevron_angle / 30) ** (0.646)) * (Re_channel_water2 ** 0.583) * (Pr_c2 ** (1/3))
    
    # Pr_l = PSI("D", "Q", 0, "P", pressureIn, wf) * PSI("V", "Q", 0, "P", pressureIn, wf) / PSI("L", "Q", 0, "P", pressureIn, wf)
    # Pho_ratio = PSI("D", "Q", 0, "P", pressureIn, wf) / PSI("D", "Q", 1, "P", pressureIn, wf)
    # pho_average = 1 / ((average_quality_in_out / PSI("D", "Q", 1, "P", pressureIn, wf)) + (1 - average_quality_in_out) / PSI("D", "Q", 0, "P", pressureIn, wf))
    # v_tp = pho_average * (average_quality_in_out * PSI("V", "Q", 1, "P", pressureIn, wf) / PSI("D", "Q", 1, "P", pressureIn, wf) + (1 - average_quality_in_out) * PSI("V", "Q", 0, "P", pressureIn, wf) / PSI("D", "Q", 0, "P", pressureIn, wf))
    # G_eq = G_channel_wf * (1 - average_quality_in_out + average_quality_in_out * (Pho_ratio ** 0.5))
    # Re_eq = G_eq * d_h / v_tp

    # Nu_tp = 12.47 * (Re_eq ** 0.33) * (Pr_l ** (1/3))
    # h_tp = Nu_tp * (PSI("L", "Q", 0, "P", pressureIn, wf)) / d_h
    # u_2 = 1 / ((1 / h_c2) + (1 / h_tp) + (area_catalogo * R_f_h) + (area_catalogo * R_f_c))
    # a_2 = ntu_2 * c_water_2 / u_2

    # Pr_h3 = PSI("V", "T", 273.15 + (PSI("T", "P", pressureIn, "Q", 1, wf) + c0.T.val - 273.15) / 2, "P", pressureIn, wf) * PSI("C", "T", 273.15 + (PSI("T", "P", pressureIn, "Q", 1, wf) - 273.15 + c0.T.val) / 2, "P", pressureIn, wf) / PSI("L", "T", 273.15 + (PSI("T", "P", pressureIn, "Q", 1, wf) - 273.15 + c0.T.val) / 2, "P", pressureIn, wf)
    # h_3 = 0.724 * (PSI("L", "T", 273.15 + (PSI("T", "P", pressureIn, "Q", 1, wf) - 273.15 + c0.T.val) / 2, "P", pressureIn, wf) / d_h) * ((chevron_angle / 30) ** (0.646)) * (Re_channel_wf3 ** 0.583) * (Pr_h3 ** (1/3))
    # Pr_c3 = PSI("V", "T", 273.15 + (T_cooling_mid2 + c21.T.val) / 2, "P", c20.p.val, coolingFluid1) * PSI("C", "T", 273.15 + (T_cooling_mid2 + c21.T.val) / 2, "P", c20.p.val, coolingFluid1) / PSI("L", "T", 273.15 + (T_cooling_mid2 + c21.T.val) / 2, "P", c20.p.val, coolingFluid1)
    # h_c3 = 0.724 * (PSI("L", "T", 273.15 + (T_cooling_mid2 + c21.T.val) / 2, "P", c20.p.val, coolingFluid1) / d_h) * ((chevron_angle / 30) ** (0.646)) * (Re_channel_water ** 0.583) * (Pr_c3 ** (1/3))
    # u_3 = 1 / ((1 / h_c3) + (1 / h_3) + (area_catalogo * R_f_h) + (area_catalogo * R_f_c))
    # a_3 = ntu_3 * c_min_3 / u_3

    # if Re_eq < 6000:
    #     f_tp = 6.947e5 * (Re_eq ** (-1.109)) * (Re_eq ** (-0.5))
    # else:
    #     f_tp = 31.21 * (Re_eq ** 0.04557) * (Re_eq ** (-0.5))

    # if Re_channel_wf3 < 2000:
    #     f_0 = 16 / Re_channel_wf3
    #     f_1 = (149 / Re_channel_wf3) + 0.9625
    # else:
    #     f_0 = (1.56 * math.log(Re_channel_wf3) - 3) ** (-2)
    #     f_1 = 9.75 / (Re_channel_wf3 ** 0.289)

    # f_h3 = ((math.cos(math.radians(90 - chevron_angle))) / (0.045 * math.tan(math.radians(90 - chevron_angle)) + 0.09 * math.sin(math.radians(90 - chevron_angle)) + f_0 / math.cos(math.radians(90 - chevron_angle))) + ((1 - math.cos(math.radians(90 - chevron_angle))) / (math.sqrt(3.8 * f_1)))) ** (-2)

    # if Re_channel_water < 2000:
    #     f_0 = 16 / Re_channel_water
    #     f_1 = (149 / Re_channel_water) + 0.9625
    # else:
    #     f_0 = (1.56 * math.log(Re_channel_water) - 3) ** (-2)
    #     f_1 = 9.75 / (Re_channel_water ** 0.289)

    # f_c = ((math.cos(math.radians(90 - chevron_angle))) / (0.045 * math.tan(math.radians(90 - chevron_angle)) + 0.09 * math.sin(math.radians(90 - chevron_angle)) + f_0 / math.cos(math.radians(90 - chevron_angle))) + ((1 - math.cos(math.radians(90 - chevron_angle))) / (math.sqrt(3.8 * f_1)))) ** (-2)

    # delta_p_f_h2 = (4 * f_tp * l_p * (G_channel_wf ** 2)) / ((PSI("D", "P", pressureIn, "Q", 1, wf) + PSI("D", "P", pressureIn, "Q", 0, wf)) * d_h)
    # delta_p_f_h3 = (4 * f_h3 * l_p * (G_channel_wf ** 2)) / (2 * PSI("D", "P", pressureIn, "T", 273.15 + (c0.T.val + PSI("T", "P", pressureIn, "Q", 1, wf) - 273.15) / 2, wf) * d_h)

    # delta_p_f_h = delta_p_f_h2 + delta_p_f_h3

    # delta_p_f_c = (4 * f_c * l_p * (G_channel_water ** 2)) / (2 * PSI("D", "P", c20.p.val, "T", 273.15 + (c20.T.val + c21.T.val) / 2, coolingFluid1) * d_h)

    # delta_p_p_c = 1.4 * (G_port_c ** 2) * 2 / (2 * PSI("D", "P", c20.p.val, "T", 273.15 + (c20.T.val + c21.T.val) / 2, mix))
    # delta_p_p_h = 1.4 * (G_port_h ** 2) * 2 / ((PSI("D", "P", pressureIn, "Q", 1, wf) + PSI("D", "P", pressureIn, "Q", 0, wf)))

    # delta_p_c = delta_p_f_c + delta_p_p_c
    # delta_p_h = delta_p_f_h + delta_p_p_h
    # delta_p_c_ratio = 1 - (delta_p_c / c20.p.val)
    # delta_p_h_ratio = 1 - (delta_p_h / c4.p.val)

    # #Dimensionamento tubo alhetado

    # # Re_outside = air_velocity * d_out * PSI("D", "T", 273.15 + T_amb, "P", 1.01325e5, "Air") / PSI("V", "T", 273.15 + T_amb, "P", 1.01325e5, "Air")

    # # if Re_outside >= 0.4 and Re_outside < 4:
    # #     c = 0.989
    # #     m_re = 0.33
    # # elif Re_outside >= 4 and Re_outside < 40:
    # #     c = 0.911
    # #     m_re = 0.385
    # # elif Re_outside >= 40 and Re_outside < 4000:
    # #     c = 0.683
    # #     m_re = 0.466
    # # elif Re_outside >= 4000 and Re_outside < 40000:
    # #     c = 0.193
    # #     m_re = 0.618
    # # else:
    # #     c = 0.027
    # #     m_re = 0.805

    # # Pr_outside = PSI("C", "T", 273.15 + T_amb, "P", 1.01325e5, "Air") * PSI("V", "T", 273.15 + T_amb, "P", 1.01325e5, "Air") / PSI("L", "T", 273.15 + T_amb, "P", 1.01325e5, "Air")
    # # Nu_outside = c * (Re_outside ** m_re) * (Pr_outside ** (1/3))
    # # h_outside = Nu_outside * PSI("L", "T", 273.15 + T_amb, "P", 1.01325e5, "Air") / d_out

    # # m = math.sqrt(2 * h_outside / (k_fin * t_fin))
    # # mr1 = m * r1_fin
    # # mr2 = m * r2c_fin
    # # I0_mr1 = iv(0, mr1)
    # # I0_mr2 = iv(0,mr2)
    # # I1_mr1 = iv(1, mr1)
    # # I1_mr2 = iv(1, mr2)
    # # K0_mr1 = kv(0, mr1)
    # # K0_mr2 = kv(0, mr2)
    # # K1_mr1 = kv(1, mr1)
    # # K1_mr2 = kv(1, mr2)

    # # Pr_wf_1 = PSI("D", "T", 273.15 + (c1.T.val + PSI("T", "Q", 1, "P", c1.p.val, wf) - 273.15) / 2, "P", c1.p.val, wf) * PSI("V", "T", 273.15 + (c1.T.val + PSI("T", "Q", 1, "P", c1.p.val, wf) - 273.15) / 2, "P", c1.p.val, wf) / PSI("L", "T", 273.15 + (c1.T.val + PSI("T", "Q", 1, "P", c1.p.val, wf) - 273.15) / 2, "P", c1.p.val, wf)
    # # Re_stage1 = 4 * c0.m.val / (math.pi * d_in * PSI("V", "T", 273.15 + (c1.T.val + PSI("T", "Q", 1, "P", c1.p.val, wf) - 273.15) / 2, "P", c1.p.val, wf))
    # # Nu_stage1 = 0.023 * (Re_stage1 ** 0.8) * (Pr_wf_1 ** 0.3)
    # # h_stage1 = Nu_stage1 * PSI("L", "T", 273.15 + (c1.T.val + PSI("T", "Q", 1, "P", c1.p.val, wf) - 273.15) / 2, "P", c1.p.val, wf) / d_in

    # # if c2.x.val ==-1:
    # #     Pr_wf_3 = PSI("D", "T", 273.15 + (c2.T.val + PSI("T", "Q", 0, "P", c1.p.val, wf) - 273.15) / 2, "P", c1.p.val, wf) * PSI("V", "T", 273.15 + (c2.T.val + PSI("T", "Q", 0, "P", c1.p.val, wf) - 273.15) / 2, "P", c1.p.val, wf) / PSI("L", "T", 273.15 + (c2.T.val + PSI("T", "Q", 0, "P", c1.p.val, wf) - 273.15) / 2, "P", c1.p.val, wf)
    # #     Re_stage3 = 4 * c0.m.val / (math.pi * d_in * PSI("V", "T", 273.15 + (c2.T.val + PSI("T", "Q", 0, "P", c1.p.val, wf) - 273.15) / 2, "P", c1.p.val, wf))
    # #     Nu_stage3 = 0.023 * (Re_stage3 ** 0.8) * (Pr_wf_3 ** 0.3)
    # #     h_stage3 = Nu_stage3 * PSI("L", "T", 273.15 + (c2.T.val + PSI("T", "Q", 0, "P", c1.p.val, wf) - 273.15) / 2, "P", c1.p.val, wf) / d_in
    # # else:
    # #     h_stage3 = 0

    # # eff_fin = (2 * r1_fin * (K1_mr1 * I1_mr2 - I1_mr1 * K1_mr2)) / (m *((r2c_fin ** 2) - (r1_fin ** 2)) * (K0_mr1 * I1_mr2 + I0_mr1 * K1_mr2))

    # # if c2.x.val == -1:
    # #     average_quality_in_out = 0.5
    # # else:
    # #     average_quality_in_out = (1 - c2.x.val) / 2
        
    # # Re_tube_l = 4 * c0.m.val * (1 - average_quality_in_out) / (math.pi * d_in * PSI("V", "T", 273.15 + c4.T.val, "Q", 0, wf))
    # # Pr_tube_l = PSI("D", "T", 273.15 + c4.T.val, "Q", 0, wf) * PSI("V", "T", 273.15 + c4.T.val, "Q", 0, wf) / PSI("L", "T", 273.15 + c4.T.val, "Q", 0, wf)
    # # Martinelli = (((1 - average_quality_in_out) / average_quality_in_out) ** 0.9) * ((PSI("D", "T", 273.15 + c4.T.val, "Q", 1, wf) / PSI("D", "T", 273.15 + c4.T.val, "Q", 0, wf)) ** 0.5) * ((PSI("V", "T", 273.15 + c4.T.val, "Q", 0, wf) / PSI("V", "T", 273.15 + c4.T.val, "Q", 1, wf)) ** 0.1)
    # # Nu_stage2 = 0.023 * (Re_tube_l ** 0.8) * (Pr_l ** 0.4) * (1 + 2.22 / (Martinelli ** 0.89))
    # # h_stage2 = Nu_stage2 * PSI("L", "T", 273.15 + c4.T.val, "Q", 0, wf) / d_in

    # # Q_stage1_L = (c1.T.val - T_amb) / ((1 / (2 * math.pi * r1_fin * h_stage1)) + (math.log(r2_fin / r1_fin) / (2 * math.pi * k_tube)) + ((S_fin + t_fin) / (h_outside * (2 * math.pi * r2_fin * S_fin + A_fin * eff_fin))))
    # # length_1 = (c1.h.val * 1000 - PSI("H", "P", c1.p.val, "Q", 1, wf)) * c0.m.val / Q_stage1_L
    # # Q_stage2_L = (PSI("T", "Q", 1, "P", c1.p.val, wf) - 273.15 - T_amb) / ((1 / (2 * math.pi * r1_fin * h_stage2)) + (math.log(r2_fin / r1_fin) / (2 * math.pi * k_tube)) + ((S_fin + t_fin) / (h_outside * (2 * math.pi * r2_fin * S_fin + A_fin * eff_fin))))
    # # length_2 = (PSI("H", "P", c1.p.val, "Q", 1, wf) - PSI("H", "P", c1.p.val, "Q", 0, wf)) * c0.m.val / Q_stage2_L
    # # if c2.x.val == -1:
    # #     Q_stage3_L = (PSI("T", "Q", 0, "P", c1.p.val, wf) - 273.15 - T_amb) / ((1 / (2 * math.pi * r1_fin * h_stage3)) + (math.log(r2_fin / r1_fin) / (2 * math.pi * k_tube)) + ((S_fin + t_fin) / (h_outside * (2 * math.pi * r2_fin * S_fin + A_fin * eff_fin))))
    # #     length_3 = ((PSI("H", "P", c1.p.val, "Q", 0, wf) - c2.h.val * 1000) * c0.m.val) / Q_stage3_L
    # # else:
    # #     length_3 = 0
    # # length_total = length_1 + length_2 + length_3

    # while abs((Cooling2.pr2.val - delta_p_h_ratio) / Cooling2.pr2.val) > 0.00001:
    #     Cooling2.pr2.val = (Cooling2.pr2.val + delta_p_h_ratio) / 2
    #     heatPump_cooling.solve(mode='design')

    # while abs((Cooling2.pr1.val - delta_p_c_ratio) / Cooling2.pr1.val) > 0.00001:
    #     Cooling2.pr1.val = (Cooling2.pr1.val + delta_p_c_ratio) / 2
    #     heatPump_cooling.solve(mode='design')

    # # d_cap = 6.35e-3
    # # a_cap = math.pi * (d_cap ** 2) / 4
    # # delta_L_2 = 0
    # # t1 = c2.T.val
    # # p1 = c2.p.val
    # # v1 = 4 * massRateWF / (math.pi * PSI("D", "P", c2.p.val, "H", 1000 * c2.h.val, wf) * (d_cap ** 2))
    # # t2 = t1 - 0.1
    # # c31.set_attr(x = None, T = t2)
    # # heatPump_cooling.solve(mode='design')
    # # x = c31.x.val
    # # p2 = c31.p.val
    # # v2 =  4 * massRateWF / (math.pi * PSI("D", "T", 273.15 + t2, "Q", c31.x.val, wf) * (d_cap ** 2))
    # # vm = (v1 + v2) / 2
    # # delta_L_inc2 = ((p1 - p2) * (a_cap / massRateWF) + (v1 - v2)) * (2 * d_cap / vm)
    # # delta_L_2 = delta_L_2 + delta_L_inc
    # # while abs(p2 - c4.p.val) / c4.p.val >= 0.01:
    # #     t1 = t2
    # #     t2 = t1 - 0.1
    # #     p1 = p2
    # #     v1 = v2
    # #     c31.set_attr(x = None, T = t2)
    # #     heatPump_cooling.solve(mode='design')
    # #     x = c31.x.val
    # #     p2 = c31.p.val
    # #     v2 =  4 * massRateWF / (math.pi * PSI("D", "T", 273.15 + t2, "Q", c31.x.val, wf) * (d_cap ** 2))
    # #     vm = (v1 + v2) / 2
    # #     delta_L_inc_2 = ((p1 - p2) * (a_cap / massRateWF) + (v1 - v2)) * (2 * d_cap / vm)
    # #     delta_L_2 = delta_L_2 + delta_L_inc_2
    
    # result_dict = {}    
    # result_dict.update({Cooling2.label: Cooling2.get_plotting_data()[2]})
    # result_dict.update({compressor.label: compressor.get_plotting_data()[1]})
    # result_dict.update({Cooling1.label: Cooling1.get_plotting_data()[1]})
    # result_dict.update({valve1.label: valve1.get_plotting_data()[1]})
    # result_dict.update({valve2.label: valve2.get_plotting_data()[1]})
    # result_dict.update({valve_three_wayA: valve_three_wayA.get_plotting_data()[1]})
    # result_dict.update({valve_three_wayB: valve_three_wayB.get_plotting_data()[1]})
    # result_dict.update({needle1.label: needle1.get_plotting_data()[1]})
    # result_dict.update({needle2.label: needle2.get_plotting_data()[1]})
    # result_dict.update({circuit_filter.label: circuit_filter.get_plotting_data()[1]})
    # result_dict.update({safety.label: safety.get_plotting_data()[1]})

    # diagram = FluidPropertyDiagram(wf)
    # diagram.set_unit_system(T='°C', p='Pa', h='kJ/kg')

    # for key, data in result_dict.items():
    #     result_dict[key]['datapoints'] = diagram.calc_individual_isoline(**data)

    # diagram.calc_isolines()

    # fig, ax = plt.subplots(1, figsize=(16, 10))
    # diagram.draw_isolines(fig, ax, 'logph', x_min= c3.h.val - 350, x_max= c1.h.val + 100, y_min= 1e4  , y_max= c1.p.val * 10)

    # for key in result_dict.keys():
    #     datapoints = result_dict[key]['datapoints']
    #     ax.plot(datapoints['h'], datapoints['p'], color='#ff0000')
    #     ax.scatter(datapoints['h'][0], datapoints['p'][0], color='#ff0000')

    # fig.savefig(str(wf) + '_logph_cooling.svg')

    # heatPump_cooling.save('hp_design_cooling')
    # heatPump_cooling.print_results()
   
    # print('O EER da bomba de calor com ' + str(wf)+ ' é ' + str(round(eer, 3)))
    # print("O rácio de compressão no arrefecimento é " + str(round(pressureRatio_cooling, 2)))
    # print("No arrefecimento, a velocidade do fluído de trabalho à saída do compressor é " + str(round(4 * massRateWF / (math.pi * PSI("D", "T", 273.15 + c1.T.val, "P", pressureOut, wf) * (0.0127 ** 2)), 2)) + " m/s")
    # print("Área necessária (permutador de placas): " + str(round(a_1 + a_2 + a_3, 3)) + " m2")
    # print("Área catalogo (permutador de placas): " + str(round(area_catalogo, 3)) + " m2")
    # print("Perda de carga no lado da água: " + str(round(delta_p_c / 1000, 4)) + " kPa")
    # print("Rácio de pressão no lado da água: " + str(round(1 - delta_p_c / c20.p.val, 3)))
    # print("Perda de carga no lado do fluído de trabalho: " + str(round(delta_p_h / 1000, 3)) + " kPa")
    # print("Rácio de pressão no lado da fluído de trabalho após permutador de placas: " + str(round(1 - delta_p_h / c4.p.val, 3)))
    # # print("O tubo alhetado tem um comprimento total de: " + str(round(length_total, 2)) + " m")

    # print("\n")

    # print("O caudal volúmico do fluído de trabalho é " + str(round(massRateWF / PSI("D", "T", 273.15 + c1.T.val, "P", pressureOut, wf) * 3600, 2)) + " m3/h") #Caudal volúmico do fluído de trabalho em kg/h
    # print("O caudal volúmico da água no aquecimento é " + water_rate_heating_v + " m3/h")
    # print("A velocidade da água no aquecimento é " + water_rate_heating + " m/s")
    # print("O caudal volúmico da água no arrefecimento é " + str(round(c20.m.val / PSI("D", "T", 273.15 + (c20.T.val + c21.T.val) / 2, "P", pressureMix, mix) * 3600, 2)) + " m3/h") #Caudal volúmico da água em kg/h
    # print("A velocidade da água no arrefecimento é " + str(round((c20.m.val / PSI("D", "T", 273.15 + (c20.T.val + c21.T.val) / 2, "P", pressureMix, mix)) / (0.25 * math.pi * (port_D_cooling ** 2)), 2)) + " m/s")

    # # print("\n")

    # # print("No aquecimento o capilar terá um comprimento de " + str(round(delta_L, 3)) + " m")
    # # print("No arrefecimento o capilar terá um comprimento de " + str(round(delta_L_2, 3)) + " m")