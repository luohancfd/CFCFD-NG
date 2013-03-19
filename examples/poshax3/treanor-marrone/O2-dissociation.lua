scheme{
   update = 'chemical kinetic ODE MC'
}

D = 59380.0
U1 = D/3.0
U2 = U1

reaction{
   'O2 + M <=> O + O + M',
   fr={'MarroneTreanor', A=1.1e25, n=-2.5, T_a=59380.0, v_name='O2', U=U1},
   chemistry_energy_coupling={{species='O2', mode='vibration', model='TreanorMarrone', U=U2}},
   ec={model='from CEA curves', iT=0}
}

