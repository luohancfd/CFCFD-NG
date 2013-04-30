-- Author: Rowan J. Gollan
-- Date: 30-Apr-2013
-- Place: The University of Queensland, Brisbane, Australia
--
-- Vibrational relaxation mechanisms identified by Park for air.
--
-- Reference:
-- Park, C. (1993)
-- Review of Chemical-Kinetic Problems of Future
-- NASA Missions, I: Earth Entries
-- Journal of Thermophysics and Heat Transfer, 7:(3), pp. 385--398
--

-- High-temperature correction is the same for all relaxation times
-- so just set it once here and assign it to label: HTC_params
HTC_params = {type='Park', sigma_dash=3e-17}

-- N2 exchanges
mechanism{
   'N2 ~~ N : V-T',
   rt={'Millikan-White:HTC', a=180, b=0.0262, HTCS=HTC_params}
}

mechanism{
   'N2 ~~ O : V-T',
   rt={'Millikan-White:HTC', a=72.4, b=0.0150, HTCS=HTC_params}
}

mechanism{
   'N2 ~~ N2 : V-T',
   rt={'Millikan-White:HTC', a=221, b=0.0290, HTCS=HTC_params}
}

mechanism{
   'N2 ~~ O2 : V-T',
   rt={'Millikan-White:HTC', a=229, b=0.0295, HTCS=HTC_params}
}

mechanism{
   'N2 ~~ NO : V-T',
   rt={'Millikan-White:HTC', a=225, b=0.0293, HTCS=HTC_params}
}

-- O2 exchanges
mechanism{
   'O2 ~~ N : V-T',
   rt={'Millikan-White:HTC', a=72.4, b=0.015, HTCS=HTC_params}
}

mechanism{
   'O2 ~~ O : V-T',
   rt={'Millikan-White:HTC', a=47.7, b=0.059, HTCS=HTC_params}
}

mechanism{
   'O2 ~~ N2 : V-T',
   rt={'Millikan-White:HTC', a=134, b=0.0295, HTCS=HTC_params}
}

mechanism{
   'O2 ~~ O2 : V-T',
   rt={'Millikan-White:HTC', a=138, b=0.0300, HTCS=HTC_params}
}

mechanism{
   'O2 ~~ NO : V-T',
   rt={'Millikan-White:HTC', a=136, b=0.0298, HTCS=HTC_params}
}

-- NO exchanges : *all* assumed same as NO-NO (which is a bit dodgy)
mechanism{
   'NO ~~ (N, O, N2, O2, NO) : V-T',
   rt={'Millikan-White:HTC', a=49.5, b=0.042, HTCS=HTC_params}
}
