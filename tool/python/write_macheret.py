# -*- coding: utf-8 -*-
"""
Created on Mon Jun 11 11:13:52 2018

@author: Han
"""
import re


def write_macheret(v_name, atom, efficiencies, A, n, Ta, fname):
    with open(fname, 'a') as f:
        for c_name, ratio in efficiencies:
            f.write('reaction{\n')
            f.write('   \'%s + %s <=> %s + %s + %s\',\n' %
                    (v_name, c_name, atom[0], atom[1], c_name))
            f.write('   fr={\'Macheret\', A=%e, n=%f, T_a=%f, v_name=\'%s\', c_name=\'%s\'},\n' % (
                    A * ratio, n, Ta, v_name, c_name))
            f.write('   chemistry_energy_coupling={{species=\'%s\', mode=\'vibration\', model=\'Macheret\',A=%e, n=%f, T_d=%f ,c_name=\'%s\'}},\n' %
                    (v_name, A * ratio, n, Ta, c_name))
            f.write('   ec={model=\'from CEA curves\',iT=0}\n}\n')


def parse_efficiencies(eff_str):
    pairs = re.findall('([\w_]+\s*=\s*[\d\.\/]+)', eff_str)
    pairs = [i.split('=') for i in pairs]
    pairs = [(i[0].replace('_plus', '+'), eval(i[1])) for i in pairs]

    return pairs


def parse_table(table, fname='temp.txt'):
    eff_str = re.findall('efficiencies={(.*?)}', table)[0]
    efficiencies = parse_efficiencies(eff_str)

    reaction = re.findall('\'(.*<=>.*)?\'', table)[0]
    pre = reaction.split('<=>')[0]
    post = reaction.split('<=>')[1]
    v_name = pre.split('+')[0]
    v_name = v_name.strip()

    post = [i.strip() for i in post.split('+')]
    atom = post[0:2]

    coeff = re.findall('fr={(.*?)}', table)[0]
    A = float(re.findall('A\s*=\s*(\d+(?:\.\d*)?(?:[eE]-?\d+)?)', coeff)[0])
    n = float(re.findall(
        'n\s*=\s*([-+]?\d+(?:\.\d*)?(?:[eE]-?\d+)?)', coeff)[0])
    Ta = float(re.findall(
        'T_a\s*=\s*([-+]?\d+(?:\.\d*)?(?:[eE]-?\d+)?)', coeff)[0])
    write_macheret(v_name, atom, efficiencies, A, n, Ta, fname)


if (__name__ == '__main__'):
    tables = ['''reaction{
   'N2 + M <=> N + N + M',
   fr={'Macheret', A=7.0e21, n=-1.60, T_a=113200.0, p_name='N2', p_mode='vibration', s_p=0.5, q_name='N2', q_mode='translation'},
   efficiencies={N2=1.0,N2_plus=1.0,O2=1.0,O2_plus=1.0,NO=1.0,NO_plus=1.0,N=30.0/7.0,N_plus=30.0/7.0,O=30.0/7.0,O_plus=30.0/7.0},
   ec={model='from CEA curves',iT=0}
}''',   '''reaction{
   'O2 + M <=> O + O + M',
   fr={'Macheret', A=2.0e21, n=-1.5, T_a=59500.0, p_name='O2', p_mode='vibration', s_p=0.5, q_name='O2', q_mode='translation'},
   efficiencies={N2=1.0,N2_plus=1.0,O2=1.0,O2_plus=1.0,NO=1.0,NO_plus=1.0,N=5.0,N_plus=5.0,O=5.0,O_plus=5.0},
   ec={model='from CEA curves',iT=0}
}''',  '''reaction{
   'NO + M <=> N + O + M',
   fr={'Park', A=5.0e15, n=0.0, T_a=75500.0, p_name='NO', p_mode='vibration', s_p=0.5, q_name='NO', q_mode='translation'},
   efficiencies={N2=1.0,N2_plus=1.0,O2=1.0,O2_plus=1.0,NO=22.0,NO_plus=1.0,N=22.0,N_plus=22.0,O=22.0,O_plus=22.0},
   ec={model='from CEA curves',iT=0}
}''']
    for table in tables:
        parse_table(table, fname='temp.txt')
