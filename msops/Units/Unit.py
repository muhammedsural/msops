class Unit():
        # Define units
        # ------------------------------------------------------------------------
        # Basic Units
        m = 1.0
        kN = 1.0
        sec = 1.0

        # Length
        mm = m / 1000.0
        cm = m / 100.0
        inch = 25.4 * mm
        ft = 12.0 * inch

        # Force
        N = kN / 1000.0
        kips = kN * 4.448221615
        lb = kips / 1.0e3

        # Stress (kN/m2 or kPa)
        Pa = N / (m ** 2)
        kPa = Pa * 1.0e3
        MPa = Pa * 1.0e6
        GPa = Pa * 1.0e9
        ksi = 6.8947573 * MPa
        psi = 1e-3 * ksi

        # Mass - Weight
        tonne = kN * sec ** 2 / m
        kg = N * sec ** 2 / m
        lb = psi*inch**2

        # Gravitational acceleration
        g = 9.81*m/sec**2

        # Time
        min = 60*sec
        hr = 60*min