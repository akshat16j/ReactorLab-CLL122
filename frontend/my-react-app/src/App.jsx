import React, { useState } from 'react';
import { BrowserRouter as Router, Routes, Route, NavLink } from 'react-router-dom';
import { Bar } from 'react-chartjs-2';
import {
  Chart as ChartJS,
  CategoryScale,
  LinearScale,
  BarElement,
  Title,
  Tooltip,
  Legend
} from 'chart.js';
// Register ChartJS components
ChartJS.register(
  CategoryScale,
  LinearScale,
  BarElement,
  Title,
  Tooltip,
  Legend
);
import './App.css';

const App = () => {
  return (
    <Router>
      <div className="app">
        <Sidebar />
        <div className="main-content">
          <Routes>
            <Route path="/" element={<Home />} />
            <Route path="/theory" element={<Theory />} />
            <Route path="/calculator" element={<Calculator />} />
          </Routes>
        </div>
      </div>
    </Router>
  );
};

const Sidebar = () => {
  return (
    <div className="sidebar">
      <div className="sidebar-header">
        <img src="/assets/reactor.png" style={{width:"99px",height:"99px"}} alt="Batch Reactor Diagram"  />
        <h1>ReactorLab</h1>
      </div>
      <nav>
        <ul>
          <li>
            <NavLink 
              to="/" 
              className={({ isActive }) => isActive ? 'active' : ''}
            ><div style={{alignItems:"center", display:"flex",}}>
              <img src="/assets/home-icon-silhouette (1).png" style={{width:"26px",height:"26px", marginRight :"15px"}} alt="Batch Reactor Diagram"  />
              Home
              </div>
            </NavLink>
          </li>
          <li>
            <NavLink 
              to="/theory" 
              className={({ isActive }) => isActive ? 'active' : ''}
            ><div style={{alignItems:"center", display:"flex",}}>
              <img src="/assets/draw (1).png" style={{width:"26px",height:"26px", marginRight :"15px"}} alt="Batch Reactor Diagram"  />
              Theory
              </div>
            </NavLink>
          </li>
          <li>
            <NavLink 
              to="/calculator" 
              className={({ isActive }) => isActive ? 'active' : ''}
            ><div style={{alignItems:"center", display:"flex",}}>
              <img src="/assets/calculator (1).png" style={{width:"26px",height:"26px", marginRight :"15px"}} alt="Batch Reactor Diagram"  />
              Calculator
              </div>
            </NavLink>
          </li>
        </ul>
      </nav>
    </div>
  );
};

const Home = () => {
  return (
    <div className="page home">
      <div className="search-bar">
        {/* ICON: search-icon.png */}
        <input type="text" placeholder="Search reactors, theories, concepts" />
      </div>
      
      
      <div className="hero">
        <h2>"From Theory to Calculation — All in One Place"</h2>
      </div>
      
      <div className="reactor-types">
        <h2>Types of Reactors:</h2>
        
        <div className="reactor-card">
          <h3>Batch Reactor</h3>
          <p>A batch reactor is a closed vessel used for carrying out chemical reactions without any continuous inflow or outflow of materials during the reaction process.</p>
          <p><strong>Key Assumptions (Ideal Batch Reactor)</strong>: Perfect Mixing: Uniform composition and temperature throughout the vessel at any time. Constant Volume: Volume remains unchanged during the reaction (unless significant gas evolution/absorption occurs). Isothermal or Adiabatic: Reactor can be operated at constant temperature (isothermal) or without heat exchange (adiabatic). Homogeneous Reaction: Reaction rate depends only on concentration (and possibly temperature), not on position or time within the reactor.</p>
        </div>
        
        <div className="reactor-card">
          <h3>Continuous Tank Stirred Reactor (CSTR)</h3>
          <img src="/assets/cstr-removebg-preview.png" style={{width:"278px",height:"264px"}} alt="Batch Reactor Diagram"  />
          <p>A Continuous Stirred Tank Reactor (CSTR) is an industrial reactor where reactants are continuously fed, mixed thoroughly, and products are continuously removed. It operates at steady state, ensuring uniform composition and temperature. Required Volume of CSTR depends on conversion of reactant and rate of reaction.</p>
          <p><strong>Key Assumptions</strong>: Perfect Mixing, Steady-State Operation, Constant Density</p>
        </div>
        
        <div className="reactor-card">
          <h3>Plug Flow Reactor (PFR)</h3>
          <img src="/assets/pfr-removebg-preview.png" style={{width:"486px",height:"70px"}} alt="Batch Reactor Diagram"  />
          <p>A Plug Flow Reactor (PFR) is a continuous-flow reactor where fluids move through a cylindrical tube in discrete "plugs" with no axial mixing, enabling precise control over reaction progression.</p>
          <p><strong>Key Characteristics</strong>: Unidirectional Flow, Gradual Concentration Gradient: Reactant concentration decreases along the reactor length, maximizing reaction rates at the inlet, Steady-State Operation: Constant flow rates and stable temperature/concentration profiles ensure consistent product quality.</p>
        </div>
        
        <div className="reactor-card">
          <h3>Packed Bed Reactor (PBR)</h3>
          <img src="/assets/pbr-removebg-preview.png" style={{width:"310px",height:"160px"}} alt="Batch Reactor Diagram"  />
          <p>A Packed Bed Reactor (PBR) is a catalytic or adsorption reactor where reactants flow through a stationary bed of solid particles (e.g., catalysts, adsorbents) to achieve chemical conversion. Widely used in gas-phase reactions, PBRs balance high efficiency with operational simplicity.</p>
          <p><strong>Key Features</strong>: Steady-State Flow: Continuous operation with minimal back mixing, Pressure Drop: Governed by the Ergun equation, critical for gas-phase systems.</p>
        </div>
      </div>
    </div>
  );
};

const Theory = () => {
  return (
    <div className="page theory">

      <div className="search-bar">
        {/* ICON: search-icon.png */}
        <input type="text" placeholder="Search reactors, theories, concepts" />
      </div>
      
      <section className="introduction">
        <h2>Introduction</h2>
        <p>This report presents the design, simulation, and analysis of chemical reactors, focusing on Batch Reactors, Plug Flow Reactors (PFR), Continuously Stirred Tank Reactors (CSTRs), Plug Bed Reactors (PBR). The study includes single and multiple reactions, energy balances, and non-isothermal effects. It also explores advanced topics like levenspiel plots and the impact of various operating conditions on conversion and selectivity.</p>
      </section>
      
      <section className="single-reactions">
        <h2>Single Reactions</h2>
        
        <div className="reactor-theory">
          <h3>1. Batch Reactor</h3>
          <p><strong>Key Assumptions (Ideal Batch Research)</strong> Perfect Mixing, Constant Volume, Isothermal or Adiabatic, Homogeneous Reaction:</p>
          <p>Governing equations:</p>
          <ul>
            <li>Material Balance: Design Equations</li>
            <li>For a reaction A {'->'} products : 
              <img src="/formulas/calculator (1).png" style={{width:"174px",height:"46px"}} alt="Batch Reactor Diagram" />
            </li>
            <li>Na: moles of reactant \( A \)
            </li>
            <li>rA: rate of reaction (mol/L.s)
            </li>
            <li>V: reactor volume (L)
            </li>
          </ul>
        </div>
        
        <div className="reactor-theory">
          <h3>2. Continuous Task Stirred Reactor (CSTR)</h3>
          <p><strong>Key Assumptions</strong> Perfect Mixing, Steady-state Operation, Constant density</p> <p>Governing equations:</p>
          <p>Design Equation:
            <img src="/formulas/image 16.png" style={{width:"270px",height:"32px"}} alt="Batch Reactor Diagram"  /></p>
          <p>For first order:
          <img src="/formulas/image 17.png" style={{width:"283px",height:"45px"}} alt="Batch Reactor Diagram"  />
          </p>
          <p>For second order:
          <img src="/formulas/image 18.png" style={{width:"294px",height:"42px"}} alt="Batch Reactor Diagram"  /></p>
          <p>Energy balance for non-isothermal and adiabatic reactions:
          <img src="/formulas/image 19.png" style={{width:"332px",height:"63px"}} alt="Batch Reactor Diagram"  /></p>
          <p>T, To : Final and initial temperature</p>
          <p>ΔHrxn​ : Heat of reaction (KJ/mol)</p>
          <p>X : Conversion of limiting reactant</p>
          <p>θi​ : Molar ratio of species i to the limiting reactant</p>
          <p>Cpi​ : Heat capacity of species i (usually in J/mol·K)</p>
          
          <p><strong>Advantages and limitations:</strong></p>
          <ul>
            <li>Pros: Excellent mixing temperature control, scalability, and continuous operation</li>
            <li>Cons: Lower conversion efficiency compared to PFRs, larger reactor size for high conversion</li>
          </ul>
        </div>
        
        <div className="reactor-theory">
          <h3>3. Plug Flow Reactor (PFR)</h3>
          <p><strong>Key Characteristics</strong> Unidirectional Flow, Gradual Concentration Gradient, Steady-State Operation</p>
          <p>Governing Equations: </p>
          <p>1. <strong>Mole Balance:</strong> <img src="/formulas/image 20.png" style={{width:"487px",height:"45px"}} alt="Batch Reactor Diagram"  /> (Design Equations)</p>
          <ul>
            <li>\( V \) - Reactor Volume
              \( FAo \) - Molar flow rate of species A entering the reactor (mol/s)
              \( X \) - Conversion of species A
              \( rA \) - Rate of reaction of A (mol/L·s)

</li>
          </ul>
          <p>2. <strong>Reaction Rate Integration</strong>
          <img src="/formulas/image 21.png" style={{width:"500px",height:"117px"}} alt="Batch Reactor Diagram"  /></p>
          <ul>
            <li>\( k \) - reaction constant
              \( CA​, CAo \) -  Concentration of species A at time t (mol/L) and initial concentration
              \( τ \) - Residence time (in s)
              \( v0 \) - volumetric flow rate</li>
          </ul>
          <p>3. <strong>Energy Balance</strong>:
            <img src="/formulas/image 22.png" style={{width:"373px",height:"47px"}} alt="Batch Reactor Diagram"  /> </p>
          <p>ΔHrxn​ : Heat of reaction (kJ/mol) — negative for exothermic reactions</p>
          <p>U : Overall heat transfer coefficient (W/m²·K)</p>
          <p>Tc​ : Coolant temperature (K or °C) </p>
          <p><img src="/formulas/image 23.png" style={{width:"274px",height:"57px"}} alt="Batch Reactor Diagram"  /></p>
          <p>ε : A small parameter (often used for perturbation or strain)</p>
          <p>
          <img src="/formulas/image 24.png" style={{width:"387px",height:"55px"}} alt="Batch Reactor Diagram"  /></p>
          <p><strong>Advantages and limitations:</strong></p>
          <ul>
            <li>Pros: high conversion efficiency and product selectivity, precise temperature control via heat exchangers, sensible design for large-scale continuous production</li>
            <li>Cons: Higher upfront costs for tubular/reactor setup, Less suitable for highly viscous fluids or reactions requiring intense mixing</li>
          </ul>
        </div>
        
        <div className="reactor-theory">
          <h3>4. Plug Bed Reactor (PRR)</h3>
          <p><strong>Key features:</strong> steady state flow, continuous operation with minimal back mixing
            Pressure Drop: governed by the ergun equation, critical for gas-phase systems.</p>
          <p><strong>Governing equations:</strong></p>
          <p>1. <strong>Mole Balance</strong>: For a reaction A {'->'} products:
          <img src="/formulas/image 25.png" style={{width:"180px",height:"73px"}} alt="Batch Reactor Diagram"  /></p>
          <p>Integrated form:
          <img src="/formulas/image 26.png" style={{width:"162px",height:"64px"}} alt="Batch Reactor Diagram"  /></p>
          <ul>
            <li>\( W \) - Catalyst weight
              \( rA' \) - Rate of disappearance of A per unit weight of catalyst (mol/kg·s)</li>
          </ul>
          <p>2. <strong>Reaction Rate</strong>:
            For a first order reaction:
            <img src="/formulas/image 27.png" style={{width:"197px",height:"37px"}} alt="Batch Reactor Diagram"  /></p>
            <p><img src="/formulas/image 28.png" style={{width:"209px",height:"33px"}} alt="Batch Reactor Diagram"  /></p>
        </div>
      </section>
    </div>
  );
};



const Calculator = () => {
  // Constants
  const R = 8.314; // J/(mol·K)

  // State for form inputs
  const [formData, setFormData] = useState({
    reactorType: '',
    CAo: '1.0',
    Ko: '0.01',
    dHrxn: '-50',
    X: '0.8',
    n: '1',
    Fo: '10',
    To: '300',
    Ea: '50',
    Cp: '75',
    pressureDropFactor: '0',
    isReversible: false,
    area: '1',
    mu: '0.001',
    eps: '0.4',
    dp: '0.01',
    P0: '101325'
  });

  const [results, setResults] = useState(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [graphData, setGraphData] = useState(null);
  const [reaction, setReaction] = useState({
    reactants: [{ coefficient: 1, species: '' }],
    products: [{ coefficient: 1, species: '' }],
    gammaA0: 1.0
  });

  // Calculation functions
  const arrhenius = (k0, Ea, T) => {
    return k0 ; // Ea converted to J/mol
  };

  const rateExpression = (CA, T, k0, Ea, n, reversible = false, CA0 = null) => {
    const k_fwd = arrhenius(k0, Ea, T);
    let r = k_fwd * Math.pow(CA, n);
    if (reversible && CA0 !== null) {
      const CA_rev = CA0 - CA;
      r -= k_fwd * Math.pow(CA_rev, n);
    }
    return r;
  };

  const calculateResults = () => {
    // Parse all inputs to numbers
    const params = {
      F0: parseFloat(formData.Fo),
      CA0: parseFloat(formData.CAo),
      T0: parseFloat(formData.To),
      k0: parseFloat(formData.Ko),
      Ea: parseFloat(formData.Ea),
      dHrxn: parseFloat(formData.dHrxn),
      Cp: parseFloat(formData.Cp),
      X: parseFloat(formData.X),
      n: parseFloat(formData.n),
      delp: parseFloat(formData.pressureDropFactor),
      isReversible: formData.isReversible,
      area: parseFloat(formData.area),
      mu: parseFloat(formData.mu),
      eps: parseFloat(formData.eps),
      dp: parseFloat(formData.dp),
      P0: parseFloat(formData.P0)
    };

    const handleReactionChange = (type, index, field, value) => {
      const updated = [...reaction[type]];
      updated[index][field] = value;
      setReaction(prev => ({ ...prev, [type]: updated }));
    };
  
    const addReactionComponent = (type) => {
      setReaction(prev => ({
        ...prev,
        [type]: [...prev[type], { coefficient: 1, species: '' }]
      }));
    };
  
    const removeReactionComponent = (type, index) => {
      setReaction(prev => ({
        ...prev,
        [type]: prev[type].filter((_, i) => i !== index)
      }));
    };

    
    try {
      let result = {};
      let profiles = {};

      if (params.X === 0) {
        // For all reactor types, volume/time should be 0 when X=0
        if (formData.reactorType === 'Batch') {
          result.time = 0;
          profiles = {
            conversion: [0],
            temperature: [params.T0],
            time: [0],
            rate: [0]
          };
        } else if (formData.reactorType === 'CSTR') {
          result.volume = 0;
          profiles = {
            conversion: [0],
            temperature: [params.T0],
            rate: [0],
            volume: [0]
          };
        } else if (formData.reactorType === 'PFR') {
          result.volume = 0;
          profiles = {
            conversion: [0],
            temperature: [params.T0],
            volume: [0],
            rate: [0]
          };
        } else if (formData.reactorType === 'PBR') {
          result.volume = 0;
          result.P_out = params.P0;
          profiles = {
            conversion: [0],
            temperature: [params.T0],
            pressure: [params.P0],
            volume: [0],
            rate: [0]
          };
        }
      } else {
        // Normal calculations when X > 0
        if (formData.reactorType === 'CSTR') {
          const T_exit = params.T0 + (-params.dHrxn * params.CA0 * params.X) / (params.F0 * params.Cp);
          const r_exit = rateExpression(
            params.CA0 * (1 - params.X),
            T_exit,
            params.k0,
            params.Ea,
            params.n,
            params.isReversible,
            params.CA0
          );
          result.volume = (params.F0 * params.X) / r_exit;
          profiles = generateCSTRProfile(params);
        } 
        else if (formData.reactorType === 'Batch') {
          const sol = solveBatch(params);
          result.time = sol.time;
          profiles = sol.profiles;
        }
        else if (formData.reactorType === 'PFR') {
          const sol = solvePFR(params);
          result.volume = sol.volume;
          profiles = sol.profiles;
        }
        else if (formData.reactorType === 'PBR') {
          const sol = solvePBR(params);
          result.volume = sol.volume;
          result.P_out = sol.P_out;
          profiles = sol.profiles;
        }
      }
  
      setResults(result);
      setGraphData(profiles);
    } catch (err) {
      setError(err.message);
    }
  };

  // Numerical integration solver (Euler's method)
  const solveIVP = (ode, initialConditions, tSpan, maxSteps = 1000, event = null) => {
    const stepSize = (tSpan[1] - tSpan[0]) / maxSteps;
    let t = tSpan[0];
    let y = [...initialConditions];
    const solution = {
      t: [t],
      y: [ [...y] ]
    };

    for (let i = 0; i < maxSteps; i++) {
      const dy = ode(t, y);
      y = y.map((val, idx) => val + dy[idx] * stepSize);
      t += stepSize;

      solution.t.push(t);
      solution.y.push([...y]);

      // Check for event
      if (event && event(t, y) <= 0) {
        break;
      }
    }

    return solution;
  };

  // Reactor-specific solvers
  const solveBatch = (params) => {
    // Improved ODE solver using Runge-Kutta 4th order
    const ode = (t, y) => {
      const [X, T] = y;
      const CA = params.CA0 * (1 - X);
      
      // Calculate reaction rate with proper reverse reaction handling
      const k_fwd = arrhenius(params.k0, params.Ea, T);
      let r = k_fwd * Math.pow(CA, params.n);
      
      if (params.isReversible) {
        const k_rev = arrhenius(params.k0_rev, params.Ea_rev, T);
        const CB = params.CA0 * (params.theta_B - (params.b/params.a) * X);
        r -= k_rev * Math.pow(CB, params.m);
      }
  
      // Correct energy balance with density (ρ) if available
      const dXdt = r / params.CA0;
      const dTdt = (params.dHrxn * r) / (params.rho * params.Cp);
  
      return [dXdt, dTdt];
    };
  
    // Adaptive step-size RK4 solver
    const adaptiveRK4 = (ode, y0, t0, tMax, maxSteps, targetX) => {
      let h = 1; // Initial step size
      let t = t0;
      let y = [...y0];
      const results = { t: [t], X: [y[0]], T: [y[1]] };
  
      for (let i = 0; i < maxSteps && t < tMax; i++) {
        let k1, k2, k3, k4;
        let yNext, error;
        
        do {
          // RK4 step
          k1 = ode(t, y);
          k2 = ode(t + h/2, y.map((yi, idx) => yi + h/2 * k1[idx]));
          k3 = ode(t + h/2, y.map((yi, idx) => yi + h/2 * k2[idx]));
          k4 = ode(t + h, y.map((yi, idx) => yi + h * k3[idx]));
          
          const y4 = y.map((yi, idx) => 
            yi + h/6 * (k1[idx] + 2*k2[idx] + 2*k3[idx] + k4[idx])
          );
          
          // Error estimation
          const y2 = y.map((yi, idx) => 
            yi + h/2 * (k1[idx] + k2[idx] + k3[idx])
          );
          
          error = Math.abs(y4[0] - y2[0]);
          if (error > 1e-4) h /= 2;
        } while (error > 1e-4 && h > 1e-6);
  
        // Update values
        t += h;
        y = y4;
        h *= 2; // Double step size for next iteration
  
        // Store results
        results.t.push(t);
        results.X.push(y[0]);
        results.T.push(y[1]);
  
        // Check for target conversion
        if (y[0] >= targetX) break;
      }
  
      // Find exact time using linear interpolation
      const lastIndex = results.X.findIndex(x => x >= params.X);
      if (lastIndex === -1) throw new Error("Target conversion not achieved");
      
      const tStart = results.t[lastIndex-1];
      const tEnd = results.t[lastIndex];
      const xStart = results.X[lastIndex-1];
      const xEnd = results.X[lastIndex];
      
      const exactTime = tStart + (params.X - xStart) * (tEnd - tStart) / (xEnd - xStart);
  
      return {
        time: exactTime,
        profiles: {
          conversion: results.X,
          temperature: results.T,
          time: results.t,
          rate: results.X.map((X, i) => {
            const CA = params.CA0 * (1 - X);
            return rateExpression(CA, results.T[i], params);
          })
        }
      };
    };
  
    // Solve with adaptive RK4
    try {
      const sol = adaptiveRK4(
        ode,
        [0, params.T0],  // Initial conditions [X, T]
        0,               // Start time
        1e6,             // Max time
        1000,            // Max steps
        params.X          // Target conversion
      );
  
      return sol;
    } catch (err) {
      throw new Error("Batch reactor solution failed: ${err.message}");
    }
  };

  const solvePFR = (params) => {
    const ode = (V, y) => {
      const [X, T] = y;
      const CA = params.CA0 * (1 - X);
      const r = rateExpression(CA, T, params.k0, params.Ea, params.n, params.isReversible, params.CA0);
      return [r / (params.F0 * params.CA0), (params.dHrxn * r) / (params.F0 * params.Cp)];
    };

    const sol = solveIVP(
      ode,
      [0, params.T0],
      [0, 1e5],
      1000,
      (V, y) => params.X - y[0] // Stop when conversion reached
    );

    return {
      volume: sol.t[sol.t.length - 1],
      profiles: {
        conversion: sol.y.map(y => y[0]),
        temperature: sol.y.map(y => y[1]),
        volume: sol.t,
        rate: sol.y.map((y, i) => rateExpression(
          params.CA0 * (1 - y[0]),
          y[1],
          params.k0,
          params.Ea,
          params.n,
          params.isReversible,
          params.CA0
        ))
      }
    };
  };

  const solvePBR = (params) => {
    const ode = (V, y) => {
      const [X, T, P] = y;
      const CA = params.CA0 * (1 - X);
      const r = rateExpression(CA, T, params.k0, params.Ea, params.n, params.isReversible, params.CA0);
      const dX = r / (params.F0 * params.CA0) * (1 - params.delp * V);
      const dT = (params.dHrxn * r) / (params.F0 * params.Cp);
      const u = params.F0 / (params.area * params.CA0);
      const dPdL = -(
        150 * params.mu * Math.pow(1 - params.eps, 2) / 
        (Math.pow(params.eps, 3) * Math.pow(params.dp, 2)) * u +
        1.75 * (1 - params.eps) / 
        (Math.pow(params.eps, 3) * params.dp) * Math.pow(u, 2)
      );
      const dP = dPdL / params.area;
      return [dX, dT, dP];
    };

    const sol = solveIVP(
      ode,
      [0, params.T0, params.P0],
      [0, 1e5],
      1000,
      (V, y) => params.X - y[0] // Stop when conversion reached
    );

    return {
      volume: sol.t[sol.t.length - 1],
      P_out: sol.y[sol.y.length - 1][2],
      profiles: {
        conversion: sol.y.map(y => y[0]),
        temperature: sol.y.map(y => y[1]),
        pressure: sol.y.map(y => y[2]),
        volume: sol.t,
        rate: sol.y.map((y, i) => rateExpression(
          params.CA0 * (1 - y[0]),
          y[1],
          params.k0,
          params.Ea,
          params.n,
          params.isReversible,
          params.CA0
        ))
      }
    };
  };

  const generateCSTRProfile = (params) => {
    // For CSTR, we'll generate a profile from 0 to final conversion
    const steps = 20;
    const conversion = [];
    const temperature = [];
    const rate = [];
    const volume = [];

    for (let i = 0; i <= steps; i++) {
      const X = (params.X * i) / steps;
      const T = params.T0 + (-params.dHrxn * params.CA0 * X) / (params.F0 * params.Cp);
      const r = rateExpression(
        params.CA0 * (1 - X),
        T,
        params.k0,
        params.Ea,
        params.n,
        params.isReversible,
        params.CA0
      );
      const V = (params.F0 * params.CA0 * X) / r;

      conversion.push(X);
      temperature.push(T);
      rate.push(r);
      volume.push(V);
    }

    return { conversion, temperature, rate, volume };
  };

  const handleChange = (e) => {
    const { name, value, type, checked } = e.target;
    setFormData(prev => ({
      ...prev,
      [name]: type === 'checkbox' ? checked : value
    }));
  };

  const handleSubmit = (e) => {
    e.preventDefault();
    setLoading(true);
    setError(null);
    
    setTimeout(() => {
      try {
        calculateResults();
      } catch (err) {
        setError(err.message);
      } finally {
        setLoading(false);
      }
    }, 100);
  };

  // Show PBR-specific fields only when PBR is selected
  const showPBRFields = formData.reactorType === 'PBR';

  // Render graphs
  // Updated renderGraphs function with centered graphs
const renderGraphs = () => {
  if (!graphData) return null;

  // Function to create a single graph
  const renderGraph = (title, xData, yData, xLabel, yLabel, color) => {
    // Calculate scales
    const xMax = Math.max(...xData);
    const yMax = Math.max(...yData);
    const yMin = Math.min(...yData);
    const yRange = yMax - yMin || 1; // Prevent division by zero

    return (
      <div className="graph-card" key={title}>
        <h4>{title}</h4>
        <div className="graph-container">
          <svg width="100%" height="300" viewBox="0 0 500 300" preserveAspectRatio="xMidYMid meet">
            {/* X axis */}
            <line x1="50" y1="250" x2="450" y2="250" stroke={color} strokeWidth="2" />
            <text x="450" y="280" fill={color} textAnchor="end">{xLabel}</text>
            
            {/* Y axis */}
            <line x1="50" y1="50" x2="50" y2="250" stroke={color} strokeWidth="2" />
            <text x="40" y="40" fill={color} textAnchor="end">{yLabel}</text>
            
            {/* Grid lines and ticks */}
            {[0, 0.2, 0.4, 0.6, 0.8, 1.0].map((frac, i) => {
              const xPos = 50 + frac * 400;
              const yPos = 250 - frac * 200;
              return (
                <g key={`grid-${i}`}>
                  {/* X grid */}
                  <line 
                    x1={xPos} y1="50" 
                    x2={xPos} y2="250" 
                    stroke="#353238" 
                    strokeWidth="1" 
                    strokeDasharray="2,2"
                  />
                  <text x={xPos} y="270" fill={color} textAnchor="middle" fontSize="10">
                    {(frac * xMax).toFixed(1)}
                  </text>
                  
                  {/* Y grid */}
                  <line 
                    x1="50" y1={yPos} 
                    x2="450" y2={yPos} 
                    stroke="#353238" 
                    strokeWidth="1" 
                    strokeDasharray="2,2"
                  />
                  <text x="30" y={yPos + 5} fill={color} textAnchor="end" fontSize="10">
                    {(yMin + frac * yRange).toFixed(1)}
                  </text>
                </g>
              );
            })}
            
            {/* Data line */}
            <polyline
              fill="none"
              stroke={color}
              strokeWidth="2"
              points={xData.map((x, i) => {
                const xPos = 50 + (x / xMax) * 400;
                const yPos = 250 - ((yData[i] - yMin) / yRange) * 200;
                return `${xPos},${yPos}`;
              }).join(' ')}
            />
          </svg>
        </div>
      </div>
    );
  };

  const xLabel = formData.reactorType === 'Batch' ? 'Time (s)' : 'Volume (L)';

  return (
    <div className="results-section">
      <h2>Calculation Results</h2>
      
      <div className="results-box">
        {/* ... keep your result items ... */}
      </div>

      <h3 className="graphs-title">Reactor Profiles</h3>
      <div className="graphs-container">
        <div className="graphs-column">
          {renderGraph(
            `Conversion vs ${formData.reactorType === 'Batch' ? 'Time' : 'Volume'}`,
            formData.reactorType === 'Batch' ? graphData.time : graphData.volume,
            graphData.conversion,
            xLabel,
            'Conversion',
            '#AC89E7'
          )}
          
          {renderGraph(
            `Temperature vs ${formData.reactorType === 'Batch' ? 'Time' : 'Volume'}`,
            formData.reactorType === 'Batch' ? graphData.time : graphData.volume,
            graphData.temperature,
            xLabel,
            'Temperature (K)',
            '#E74C3C'
          )}
        </div>
        
        <div className="graphs-column">
          {renderGraph(
            "Rate of Reaction vs Conversion",
            graphData.conversion,
            graphData.rate,
            'Conversion',
            'Rate (mol/L·s)',
            '#3498DB'
          )}
          
          {formData.reactorType === 'PBR' && graphData.pressure && renderGraph(
            "Pressure vs Volume",
            graphData.volume,
            graphData.pressure,
            'Volume (L)',
            'Pressure (Pa)',
            '#2ECC71'
          )}
        </div>
      </div>
    </div>
  );
};

  return (
    <div className="calculator-page">
      <div className="calculator-container">
        <h1>ReactorLab Calculator</h1>
        
        {error && <div className="error-message">{error}</div>}
        
        <form onSubmit={handleSubmit} className="calculator-form">
          <div className="form-columns">
            {/* Left Column */}
            <div className="form-column">
              <div className="form-group">
                <label>Reactor Type</label>
                <select
                  name="reactorType"
                  value={formData.reactorType}
                  onChange={handleChange}
                  required
                >
                  <option value="">Select Reactor Type</option>
                  <option value="Batch">1. Batch Reactor</option>
                  <option value="CSTR">2. CSTR</option>
                  <option value="PFR">3. PFR</option>
                  <option value="PBR">4. PBR</option>
                </select>
              </div>
                
              <div className="form-group">
                <label>Initial concentration (CAo, mol/L)</label>
                <input
                  type="number"
                  name="CAo"
                  value={formData.CAo}
                  onChange={handleChange}
                  step="0.01"
                  min="0"
                  required
                />
              </div>

              <div className="form-group">
                <label>Rate Constant (K)</label>
                <input
                  type="number"
                  name="Ko"
                  value={formData.Ko}
                  onChange={handleChange}
                  step="0.0001"
                  min="0"
                  required
                />
              </div>

              <div className="form-group">
                <label>Enthalpy of reaction (ΔHrxn, kJ/mol)</label>
                <input
                  type="number"
                  name="dHrxn"
                  value={formData.dHrxn}
                  onChange={handleChange}
                  step="0.1"
                  required
                />
              </div>

              <div className="form-group">
                <label>Target conversion (X)</label>
                <input
                  type="number"
                  name="X"
                  value={formData.X}
                  onChange={handleChange}
                  step="0.01"
                  min="0"
                  max="1"
                  required
                />
              </div>

              <div className="form-group">
                <label>Reaction Order (n)</label>
                <input
                  type="number"
                  name="n"
                  value={formData.n}
                  onChange={handleChange}
                  step="0.1"
                  min="0"
                  required
                />
              </div>
            </div>

            {/* Right Column */}
            <div className="form-column">
              <div className="form-group">
                <label>Feed flow rate (Fo, mol/s)</label>
                <input
                  type="number"
                  name="Fo"
                  value={formData.Fo}
                  onChange={handleChange}
                  step="0.01"
                  min="0"
                  required={formData.reactorType !== 'Batch'}
                />
              </div>

              <div className="form-group">
                <label>Initial temperature (To, K)</label>
                <input
                  type="number"
                  name="To"
                  value={formData.To}
                  onChange={handleChange}
                  step="1"
                  min="0"
                  required
                />
              </div>

              <div className="form-group">
                <label>Activation energy (Ea, kJ/mol)</label>
                <input
                  type="number"
                  name="Ea"
                  value={formData.Ea}
                  onChange={handleChange}
                  step="0.1"
                  min="0"
                  required
                />
              </div>

              <div className="form-group">
                <label>Heat capacity (Cp, J/mol·K)</label>
                <input
                  type="number"
                  name="Cp"
                  value={formData.Cp}
                  onChange={handleChange}
                  step="0.1"
                  min="0"
                  required
                />
              </div>

              <div className="form-group">
                <label>Pressure drop factor (α, optional)</label>
                <input
                  type="number"
                  name="pressureDropFactor"
                  value={formData.pressureDropFactor}
                  onChange={handleChange}
                  step="0.0001"
                  min="0"
                />
              </div>

              <div className="form-group checkbox-group">
                <label>
                  <input
                    type="checkbox"
                    name="isReversible"
                    checked={formData.isReversible}
                    onChange={handleChange}
                  />
                  Is Reversible?
                </label>
              </div>

              {/* PBR Specific Fields */}
              {showPBRFields && (
                <>
                  <div className="form-group">
                    <label>Cross-sectional area (m²)</label>
                    <input
                      type="number"
                      name="area"
                      value={formData.area}
                      onChange={handleChange}
                      step="0.01"
                      min="0.01"
                      required
                    />
                  </div>

                  <div className="form-group">
                    <label>Viscosity (μ, Pa·s)</label>
                    <input
                      type="number"
                      name="mu"
                      value={formData.mu}
                      onChange={handleChange}
                      step="0.0001"
                      min="0"
                      required
                    />
                  </div>

                  <div className="form-group">
                    <label>Bed porosity (ε)</label>
                    <input
                      type="number"
                      name="eps"
                      value={formData.eps}
                      onChange={handleChange}
                      step="0.01"
                      min="0.1"
                      max="0.9"
                      required
                    />
                  </div>

                  <div className="form-group">
                    <label>Particle diameter (dp, m)</label>
                    <input
                      type="number"
                      name="dp"
                      value={formData.dp}
                      onChange={handleChange}
                      step="0.001"
                      min="0.001"
                      required
                    />
                  </div>

                  <div className="form-group">
                    <label>Initial pressure (P0, Pa)</label>
                    <input
                      type="number"
                      name="P0"
                      value={formData.P0}
                      onChange={handleChange}
                      step="1000"
                      min="0"
                      required
                    />
                  </div>
                </>
              )}
            </div>
          </div>

          <div className="form-actions">
            <button 
              type="submit" 
              className="calculate-button"
              disabled={loading}
            >
              {loading ? 'Calculating...' : 'Calculate'}
            </button>
          </div>
        </form>

        {results && (
          <div className="results-section">
            <h2>Calculation Results</h2>
            
            <div className="results-box">
              {results.volume !== undefined && (  // Check for undefined, not truthy
                <div className="result-item">
                  <span>Reactor Volume</span>
                  <span>
                    {results.volume > 0 ? results.volume.toFixed(2) : "0"} L
                  </span>
                </div>
              )}
              {results.time !== undefined && (
                <div className="result-item">
                  <span>Reaction Time</span>
                  <span>
                    {results.time > 0 ? results.time.toFixed(2) : "0"} s
                  </span>
                </div>
              )}
              {results.P_out !== undefined && (
                <div className="result-item">
                  <span>Outlet Pressure</span>
                  <span>{results.P_out.toFixed(0)} Pa</span>
                </div>
              )}
            </div>

            {renderGraphs()}
          </div>
        )}
      </div>
    </div>
  );
};

export default App;