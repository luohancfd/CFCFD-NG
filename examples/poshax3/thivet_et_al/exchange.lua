-- VT exchanges --
mechanism{
   'O2 ~~ (O2, N2) : V-T',
   rt={'SSH-VT'}
}

mechanism{
   'N2 ~~ (O2, N2) : V-T',
   rt={'SSH-VT'}
}

-- VV exchanges --
mechanism{
   'O2 ~~ N2 : V-V THO',
   rt={'SSH-VV'}
}

mechanism{
   'N2 ~~ O2 : V-V THO',
   rt={'SSH-VV'}
}


