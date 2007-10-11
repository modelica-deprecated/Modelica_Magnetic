package Modelica_Magnetic "Library for modelling of electromagnetic devices with lumped magnetic networks"
package UsersGuide "Users Guide"

  annotation (DocumentationClass=true, Documentation(info="<HTML>
<p>
This library contains components for modelling of electromagnetic devices with lumped magnetic networks. Those models are suited for both rough design of the magnetic subsystem of a device as well as for efficient dynamic simulation at system level together with neighbouring subsystems. At present, components and examples for modelling of <i>translatory</i> electromagnetic and electrodynamic actuators are provided. If needed, these components can be adapted to network modellling of <i>rotational</i> electrical machines.
</p>
<p>
This users guide gives a short introduction to the underlying concept of <b>magnetic flux tubes</b>, summarizes the calculation of magnetic <b>reluctance forces</b> from lumped magnetic network models and lists <b>reference literature</b>.
</p>
<p>
<a href=\"Modelica_Magnetic.Examples\">Examples</a> illustrates the usage of magnetic network models with simple models from different fields of application.
</HTML>"));

  class FluxTubeConcept "Flux tube concept"

    annotation (Documentation(info="<HTML>
<h3><font color=\"#008000\" size=5>Overview of the Conept of Magnetic Flux Tubes</font></h3>
<p>
Following below, the concept of magnetic flux tubes is outlined in short. For a detailed description of flux tube elements, please have a look at the listed literature. Magnetic flux tubes enable for modeling of magnetic fields with lumped networks. The figure below and the following equations illustrate the transition from the original magnetic field quantities described by <i>Maxwell</i>'s equations to network elements with a flow variable and an across variable:
<br>
<dl>
<dd>
<p><IMG SRC=\"../Images/Magnetic/UsersGuide/FluxTubeConcept/magnetic_flux_tube_schematic.png\" ALT=\"Magnetic flux tube\"></p>
</dd>
</dl>
<br>
For a region with an approximately homogeneous distribution of the magnetic field strength <b>H</b> and the magnetic flux density <b>B</b> through cross sectional area <i>A</i> at each length coordinate <i>s</i>  (<i>A</i> perpendicular to the direction of the magnetic field lines), a magnetic reluctance <i>R<sub>m</sub></i> can be defined:<br>
<dl>
<dd>
<p><IMG SRC=\"../Images/Magnetic/UsersGuide/FluxTubeConcept/eq_transition_reluctance_flowAcross_IntegralQuantities.png\" ALT=\"Transition from field quantities to flow- and across variables\"></p>
</dd>
</dl>
With the definition of the magnetic potential difference <i>V<sub>mag</sub></i> as an across variable and the magnetic flux <i>&Phi;</i> as flow variable, a reluctance element <i>R<sub>m</sub></i> can be defined similar to resistive network elements in other physical domains. Using <i>Maxwell</i>'s constitutive equation<br>
<dl>
<dd>
<p><IMG SRC=\"../Images/Magnetic/UsersGuide/FluxTubeConcept/eq_MaxwellConstitutive.png\" ALT=\"Maxwell's constitutive equation\"></p>
</dd>
</dl>
the general formula for the calculation of a magnetic reluctance <i>R<sub>m</sub></i> from its geometric and material properties is:<br>
<dl>
<dd>
<p><IMG SRC=\"../Images/Magnetic/UsersGuide/FluxTubeConcept/eq_reluctance_general.png\" ALT=\"General formula for calculation of a magnetic reluctance\"></p>
</dd>
</dl>
For a prismatic or cylindrical volume of length <i>l</i> and cross sectional area <i>A</i> with the magnetic flux entering and leaving the region through its end planes, the above equation simplifies to:<br>
<dl>
<dd>
<p><IMG SRC=\"../Images/Magnetic/UsersGuide/FluxTubeConcept/eq_reluctance_prismatic.png\" ALT=\"Magnetic reluctance of a prismatic or cylindrical volume\"></p>
</dd>
</dl>
<p>
Similar equations can be derived for other geometries. In cases where a direct integration is not possible, the reluctance can be calclulated on base of average length, average cross sectional area and volume <i>V</i> respectively:
</p>
<dl>
<dd>
<p><IMG SRC=\"../Images/Magnetic/UsersGuide/FluxTubeConcept/eq_reluctanceFromAverageGeometry.png\" ALT=\"Reluctance calculation from average geometric quantities\"></p>
</dd>
</dl>
<p>
Network elements for sources of a magnetic potential difference or magnetomotive force, i.e. coils or permanent magnets can be formulated as well. The resulting magnetic network models of actuators reflect the main dimensions of these devices as well as the normally nonlinear characteristics of their magnetically active materials.
</p>

</HTML>
"));
  equation

  end FluxTubeConcept;

  class ReluctanceForceCalculation "Reluctance forces"

    annotation (Documentation(info="<HTML>
<h3><font color=\"#008000\" size=5>Calculation of reluctance forces from lumped magnetic network models</font></h3>

<p>
Generally, the thrust <i>F</i> developed by a translatory electro-magneto-mechanical actuator (similar for the rotational case with torque and angular position) is equal to the change of magnetic co-energy <i>W<sub>m</sub><sup>*</sup></i> with armature position <i>x</i> according to
<dl>
<dd>
<p><IMG SRC=\"../Images/Magnetic/UsersGuide/ReluctanceForceCalculation/eq_CoEnergy_general.png\" ALT=\"Equation for force calculation from change of magnetic co-energy with armature position\"></p>
</dd>
</dl>
(<i>&Psi;</i> flux linkage, <i>i</i> actuator current). In lumped magnetic network models, the above equation simplifies to
<dl>
<dd>
<p><IMG SRC=\"../Images/Magnetic/UsersGuide/ReluctanceForceCalculation/eq_forceFromPermeance_network.png\" ALT=\"Equation for force calculation in lumped magnetic network models\"></p>
</dd>
</dl>
where <i>n<sub>linear</sub></i> is the number of flux tube elements with constant relative permeability that change its permeance <i>G<sub>m i</sub></i> with armature position (index <i>i</i>), <i>V<sub>mag i</sub></i> the magnetic voltage across each respective flux tube and <i>dG<sub>m i</sub>/dx</i> the derivative of the respective permeances with respect to armature position. Transition from the general formula based on magnetic co-energy to the latter one is outlined in <a href=\"Literature\">[3]</a> for the reciprocal of the permeance, i.e. for the magnetic reluctance <i>R<sub>m</sub></i>. Note that
<dl>
<dd>
<p><IMG SRC=\"../Images/Magnetic/UsersGuide/ReluctanceForceCalculation/eq_transition_forceReluctancePermeance.png\" ALT=\"Transition from force calculation based on reluctance to calculation based on permeance\"></p>
</dd>
</dl>
with <i>&Phi;<sub>mag i</sub></i> being the magnetic flux through each respective flux tube element.
</p>
<p>
Flux tube elements with <i>nonlinear</i> material characteristics <i>&mu;<sub>r</sub></i>(<i>B</i>) in magnetic network models do not restrict the usability of the above equation. However, it is required that these nonlinear flux tube elements do not change its shape with armature motion (e.g. portion of a solenoid plunger where the magnetic flux passes through in axial direction). This limitation is not a strong one, since the permeance of nonlinear, but highly permeable ferromagnetic flux tube elements and its change with armature position compared to that of air gap flux tubes can be neglected in most cases. Because of this constraint, the dimensions of possibly nonlinear flux tube elements in sub-package <a href=\"Modelica_Magnetic.FluxTube.FixedShape\">FluxTube.FixedShape</a> are fixed, whereas the dimension in direction of motion of the linear flux tube elements in sub-package <a href=\"Modelica_Magnetic.FluxTube.Force\">FluxTube.Force</a> can vary during simulation. For the flux tubes defined in this package with their rather simple shapes, the derivative <i>dG<sub>m</sub>/dx</i> is given analytically. For more complex shapes and variations of dimensions with armature motion, it must be provided analytically during model development, preferably by extending the partial model <a href=\"Modelica_Magnetic.FluxTube.Force.PartialForce\">FluxTube.Force.PartialForce</a>.<br>
</p>
<p>
The sub-package <a href=\"Modelica_Magnetic.FluxTube.Leakage\">FluxTube.Leakage</a> contains flux tube shapes typical for leakage flux around prismatic or cylindrical poles. Since the permeance of these flux tubes does not change with armature position, they do not contribute to a reluctance actuator's thrust.
</p>


</HTML>
"));
  equation

  end ReluctanceForceCalculation;

  class Literature "Literature"

    annotation (Documentation(info="<html>
<h3><font color=\"#008000\" size=5>Literature</font></h3>
<ul>
<li>
A first realisation of the Modelica Magnetic library is described in:
<dl>
<dt>[1] B&ouml;drich, Th.; Roschke, Th.:</dt>
<dd> <b>A Magnetic Library for Modelica</b>.
     Modelica 2005 Conference, Hamburg, Germany,
     pp. 559-565, March 7-8, 2005.
     Download from:
     <a href=\"http://www.modelica.org/events/Conference2005/online_proceedings/Session7/Session7a3.pdf\">http://www.modelica.org/events/Conference2005/online_proceedings/Session7/Session7a3.pdf</a>
     </dd>
</dl>
</li>
<li>
The method of magnetic flux tubes as well as derivation of the permeance of many flux tube shapes is explained in detail in:
<dl>
<dt>[2] Roters, H.:</dt>
<dd> <b>Electromagnetic Devices</b>.
New York: John Wiley & Sons 1941
</dd>
</dl>
</li>
<li> Structure, properties, applications and design of translatory electromagnetic (i.e. reluctance type) actuators are thoroughly described in:
<dl>
<dt>[3] Kallenbach, E.; Eick, R.; Quendt, P.; Str&ouml;hla, T.; Feindt, K.; Kallenbach, M.:</dt>
<dd><b>Elektromagnete: Grundlagen, Berechnung, Entwurf und Anwendung</b>.
2nd ed., Wiesbaden: B.G. Teubner 2003
<br>&nbsp;</dd>
<dt>[4] Roschke, T.:</dt
<dd><b>Entwurf geregelter elektromagnetischer Antriebe f&uuml;r Luftsch&uuml;tze</b>.
    Fortschritt-Berichte VDI, Reihe 21, Nr. 293, D&uuml;sseldorf: VDI-Verlag 2000</dd>
</dl>
</li>
<li>
Application of the method of magnetic flux tubes to the design of rotational electrical machines is explained for example in:
<dl>
<dt>[5] Hendershot, J.R. Jr.; Miller, T.J.E.:</dt>
<dd> <b>Design of Brushless Permanent-Magnet Motors</b>.
Magna Physics Publishing and Oxford University Press 1994
</dd>
</dl>
</li>

</ul>

</html>
"));
  end Literature;

  class Contact "Contact"

    annotation (Documentation(info="<html>
<h3><font color=\"#008000\" size=5>Contact</font></h3>

<dl>
<dt>
<b>Main Author:</b></dt>
<br>
<dd>
    <a href=\"http://www.ifte.de/mitarbeiter/boedrich.html\">Thomas B&ouml;drich</a><br>
    Dresden University of Technology<br>
    Institute of Electromechanical and Electronic Design<br>
    01062 Dresden, Germany<br>
    Phone: ++49 - 351 - 463 36296<br>
    Fax: ++49 - 351 - 463 37183<br>
    email: <A HREF=\"mailto:Thomas.Boedrich@mailbox.tu-dresden.de\">Thomas.Boedrich@mailbox.tu-dresden.de</A><br></dd>
</dl>
<br>
<p><b>Acknowledgements:</b></p>
<ul>

<li>
The magnetisation characteristics of the included soft magnetic materials were compiled and measured respectively by Thomas Roschke, now with Saia-Burgess Dresden GmbH, Dresden, Germany. Provision of this data is highly appreciated. He also formulated the approximation function used for description of the magnetisation characteristics of these materials.
</li>

<li>
André Klick of then Dresden University of Technology, Dresden, Germany gave valuable support on the implementation of this library. His contribution is highly appreciated, too.
</li>

</ul>

</html>
"));
  end Contact;

end UsersGuide;


  package Examples
  "Illustration of component usage with simple models of various devices"

    model SaturatedInductor
    "Inductor with saturation in the ferromagnetic core"

      SI.MagneticFieldStrength H_fe
      "Magnetic field strength inside core for plot of characteristics B(H)";
      Sources.ElectroMagneticConverter coil(w=1000, c_coupl=0.85)
      "Inductor coil"
        annotation (extent=[10,0; 30,20]);
      MagneticGround magGround      annotation (extent=[40,-30; 60,-10]);
      Modelica_Magnetic.FluxTube.FixedShape.Cuboid R_mFe(
        l=4*0.1,
        a=0.04,
        b=0.04,
        redeclare record Material =
          Modelica_Magnetic.Material.SoftMagnetic.ElectricSheet.M350_50A)
      "Reluctance of ferromagnetic inductor core"
        annotation (extent=[40,0; 60,20], rotation=270);

      annotation (
          experiment(StopTime=0.1, Tolerance=1e-007),
        Documentation(info="<html>
<p>
This model demonstrates the effects of non-linear magnetisation characteristics of soft magnetic materials (hysteresis neglected). A sinusoidal voltage is applied to an inductor with a closed ferromagnetic core of rectangular shape. Set the <b>tolerance</b> to <b>1e-7</b>, <b>simulate for 0.1 s</b> and plot for example:
</p>
<pre>
    coil.i vs. time           // non-harmonic current due to saturation of the core material
    R_mFe.my_r vs. R_mFe.B    // relative permeability vs. flux density inside core
    R_mFe.B vs. H_fe          // magnetisation curve B(H); hysteresis neglected
</pre>
<p>
The magnetisation characteristics of the flux tube element representing the ferromagnetic core can easily be changed from simplified linear behaviour (nonLinearPermeability set to false and R_mFe.my_rConst set to a positive value, preferably my_rConst >> 1) to non-linear behaviour (e.g. selection of one of the electric sheets in <a href=\"Material.SoftMagnetic\">Material.SoftMagnetic</a> with nonLinearPermeability set to true). This enables for convenient inital design of magnetic circuits with linear material characteristics prior to simulation with non-linear behaviour. </p>
<p>
N.B.:<br>
If the supply voltage has a zero-crossing when applied to the inductor at time t=0 (i.e. u_source.phase set to zero instead of pi/2), then the inrush current that is typical for switching of inductive loads can be observed.
</p>
</html>"),Coordsys(
          extent=[-100, -100; 100, 100],
          grid=[2, 2],
          component=[20, 20]),
          Icon(
          Rectangle(extent=[-100, -100; 80, 50], style(fillColor=7)),
          Polygon(points=[-100, 50; -80, 70; 100, 70; 80, 50; -100, 50], style(
                fillColor=7)),
          Polygon(points=[100, 70; 100, -80; 80, -100; 80, 50; 100, 70], style(
                fillColor=7)),
          Text(
            extent=[-100,70; 100,-130],
            style(color=3, rgbcolor={0,0,255}),
            string="1"),
          Text(
            extent=[-120,130; 120,70],
            string="%name",
            style(color=3))),
        Diagram);

      Modelica.Electrical.Analog.Basic.Ground ground
        annotation (extent=[-70,-30; -50,-10]);
      Modelica.Electrical.Analog.Basic.Resistor R_coil(R=10)
      "Inductor coil's resistance"      annotation (extent=[-19,10; 1,30]);
      Modelica.Electrical.Analog.Sources.SineVoltage u_source(
        freqHz=50,
        V=1000,
        phase=pi/2) "Voltage applied to inductor"
                annotation (extent=[-70,0; -50,20], rotation=270);
    equation
      connect(coil.p_mag, R_mFe.p)          annotation (points=[30,16; 30,20;
            50,20],  style(color=45, rgbcolor={255,127,0}));
      connect(coil.n_mag, magGround.p)                annotation (points=[30,4; 30,
            0; 50,0; 50,-10],style(color=45, rgbcolor={255,127,0}));
      connect(u_source.p, R_coil.p) annotation (points=[-60,20; -19,20], style(
            color=3, rgbcolor={0,0,255}));
      connect(R_coil.n, coil.p_el)           annotation (points=[1,20; 10,20;
            10,16],  style(color=3, rgbcolor={0,0,255}));
      connect(coil.n_el, u_source.n)           annotation (points=[10,4; 10,0;
            -60,0], style(color=3, rgbcolor={0,0,255}));
      connect(u_source.n, ground.p)      annotation (points=[-60,0; -60,-10],
          style(color=3, rgbcolor={0,0,255}));
      connect(magGround.p, R_mFe.n)     annotation (points=[50,-10; 50,0],
                   style(color=45, rgbcolor={255,127,0}));
      H_fe = R_mFe.B/(my_0 * R_mFe.my_r);
    end SaturatedInductor;
    extends Modelica.Icons.Library;
    package ElectrodynamicActuator
    "Two translatory electrodynamic actuator models of different modelling depth and their usage"

      model MotorConstantModel
      "Simple behavioural actuator model for system simulation"
        extends Modelica_Magnetic.Interfaces.ElectromechanicalActuator;

        parameter Real c_m = 1 "Motor constant in N/A and V/(m/s) respectively";
        parameter SI.Mass m_a = 0.1 "Armature mass";
        parameter SI.Position x_min = 0 "Minimum armature position";
        parameter SI.Position x_max = 0.1 "Maximum armature position";
        parameter SI.Resistance R_coil = 1 "Coil resistance";
        parameter SI.Inductance L_coil = 0.1 "Coil inductance";

        SI.Position x(start=x_min, stateSelect=StateSelect.prefer)
        "Armature position, alias for flange position";

        Modelica.Electrical.Analog.Sensors.CurrentSensor iSensor "coil current"
          annotation (extent=[-80,10; -60,-10],  rotation=90);
        Modelica.Electrical.Analog.Sources.SignalVoltage u_i
        "Induced back-emf due to armature motion"
          annotation (extent=[-80,-40; -60,-20], rotation=270);
        Modelica.Mechanics.Translational.Force MotorForce
          annotation (extent=[10,-10; 30,10]);
        Modelica.Electrical.Analog.Basic.Resistor R(R=R_coil) "Coil resistance"
          annotation (extent=[-90,50; -70,70]);
        Modelica.Mechanics.Translational.Sensors.SpeedSensor vSensor
        "armature velocity"
          annotation (extent=[10,-40; 30,-20],
                                             rotation=180);
        Modelica_Magnetic.Utilities.TranslatoryArmature armature(
          m=m_a,
          x_max=x_max,
          x_min=x_min) "Armature inertia with stoppers at end of stroke range"
                       annotation (extent=[60,-10; 80,10]);

        Modelica.Electrical.Analog.Basic.Inductor L(L=L_coil) "Coil inductance"
          annotation (extent=[-80,20; -60,40], rotation=270);
        Modelica.Blocks.Math.Gain c_force(k=c_m) annotation (extent=[-40,-10; -20,10]);
        Modelica.Blocks.Math.Gain c_velocity(k=c_m)
          annotation (extent=[-20,-40; -40,-20]);
      equation
        flange.s = x;

        connect(iSensor.n, u_i.p)   annotation (points=[-70,-10; -70,-20], style(
              color=3, rgbcolor={0,0,255}));
        annotation (Diagram, Icon(
            Text(
              extent=[-138,150; 74,-82],
              string="c",
              style(color=45, rgbcolor={255,128,0})),
            Text(
              extent=[-6,16; 96,-94],
              style(color=45, rgbcolor={255,128,0}),
              string="m")),
          Documentation(info="<html>
<p>
This is a simple behavioural model of a translatory electrodynamic actuator (either moving coil or moving magnet type). The electro-mechanical conversion process is described with the motor constant c_m. The model is very similar to the well-known behavioural model of a rotational DC-Machine, except that it is for translatory motion. <br>
<br>
The motor constant c_m as well as coil resistance R and inductance L are assumed to be known, e.g. from measurements or from catalogue data. Hence this model is well-suited for system simulation together with neighbouring subsystems, but not for actuator design, where the motor constant is to be found on base of the magnetic circuit's geometry, material properties and winding data. See <a href=\"MagneticCircuitModel\">MagneticCircuitModel</a> for an actuator model intended for rough design of an electrodynamic actuator. Due to identical connectors, both models can be used in system simulation, e.g. to simulate a stroke as demonstrated in <a href=\"ArmatureStroke\">ArmatureStroke</a>.
</p>
</html>"));
        connect(R.p, p) annotation (points=[-90,60; -100,60], style(
            color=3,
            rgbcolor={0,0,255},
            fillPattern=1));
        connect(u_i.n, n)   annotation (points=[-70,-40; -70,-60; -100,-60],
            style(
            color=3,
            rgbcolor={0,0,255},
            pattern=0));
        connect(armature.flange_a, vSensor.flange_a) annotation (points=[60,0; 50,0;
              50,-30; 30,-30],
                     style(
            color=58,
            rgbcolor={0,127,0},
            pattern=0,
            fillColor=47,
            rgbfillColor={255,170,85},
            fillPattern=1));
        connect(armature.flange_b, flange) annotation (points=[80,0; 100,0],
            style(
            color=58,
            rgbcolor={0,127,0},
            pattern=0,
            fillColor=47,
            rgbfillColor={255,170,85},
            fillPattern=1));
        connect(L.p, R.n)
          annotation (points=[-70,40; -70,60], style(color=3, rgbcolor={0,0,255}));
        connect(L.n, iSensor.p)
          annotation (points=[-70,20; -70,10], style(color=3, rgbcolor={0,0,255}));
        connect(c_force.y, MotorForce.f)
          annotation (points=[-19,0; 8,0], style(color=74, rgbcolor={0,0,127}));
        connect(c_velocity.y, u_i.v) annotation (points=[-41,-30; -63,-30], style(
              color=74, rgbcolor={0,0,127}));
        connect(c_velocity.u, vSensor.v)
          annotation (points=[-18,-30; 9,-30], style(color=74, rgbcolor={0,0,127}));
        connect(MotorForce.flange_b, armature.flange_a)
          annotation (points=[30,0; 60,0], style(color=58, rgbcolor={0,127,0}));
        connect(c_force.u, iSensor.i) annotation (points=[-42,0; -51,0; -51,
            6.12303e-016; -60,6.12303e-016],   style(color=74, rgbcolor={0,0,127}));
      end MotorConstantModel;

      annotation (Documentation(info="<html>
<p>
Similar to rotational DC-Motors, the electro-mechanical energy conversion of translatory electrodynamic actuators can be described with the following two converter equations:
<pre>
      F = c_m * i
    V_i = c_m * v
</pre>
with electrodynamic or <i>Lorentz</i> force F, motor constant c_m, current i, induced back-emf V_i and armature velocity v. For a moving coil actuator with a coil inside an air gap with flux density B and a total wire length l inside the magnetic field, the motor constant becomes
<pre>
    c_m = B * l
</pre>
This motor constant c_m can roughly be determined with a lumped magnetic network model of the actuator's permanent magnetic excitation system. Coil resistance R and inductance L can be calculated from geometry and material data of the magnetic circuit and the winding as well.
</p>
</html>"),     Icon(
          Rectangle(extent=[-100,-100; 80,50],   style(fillColor=7)),
          Polygon(points=[-100,50; -80,70; 100,70; 80,50; -100,50],      style(
                fillColor=7)),
          Polygon(points=[100,70; 100,-80; 80,-100; 80,50; 100,70],      style(
                fillColor=7)),
          Text(
            extent=[-100,70; 100,-130],
            style(color=3, rgbcolor={0,0,255}),
            string="2"),
          Text(
            extent=[-120,130; 120,70],
            string="%name",
            style(color=3))));

      model MagneticCircuitModel
      "Detailed actuator model for rough design of actuator and system simulation"

        extends Modelica_Magnetic.Interfaces.ElectromechanicalActuator;

        parameter SI.Radius r_core = 18.5e-3
        "Radius of ferromagnetic stator core";

        parameter SI.Thickness l_PM = 3.2e-3
        "Radial thickness of permanent magnet ring";
        parameter SI.Length t = 0.02
        "Axial length of permanent magnet ring and air gap respectively";

        parameter SI.Length t_add = 0.01
        "Additional axial clearance between permanent magnet and stator bottom side";

        parameter SI.Thickness t_bot = 0.01
        "Axial thickness of stator bottom side";

        parameter SI.Thickness l_air = 3e-3
        "Total radial length of armature air gap";
        parameter SI.Thickness l_clear = 0.5e-3
        "Radial clearance armature <=> stator (on either side of armature)";
        parameter SI.Thickness l_car = 0.5e-3
        "Radial thickness of coil carrier";

        parameter SI.Thickness t_out = 6.3e-3
        "Radial thickness of outer section of stator";

        parameter SI.Breadth w_w = 12e-3
        "Width of armature winding (axial direction)";

        parameter SI.Position x_min = w_w/2 "Minimum armature position";
        parameter SI.Position x_max = t - w_w/2 "Maximum armature position";

        parameter SI.Mass m_a = 0.08 "Armature mass";

        SI.Position x(start = x_min, stateSelect=StateSelect.prefer)
        "Armature position, alias for flange position";

        SI.Reluctance R_mTot "Estimate for total reluctance as seen by coil";
        // R_mTot is different from total reluctance as seen by permanent magnet's magnetomotive force, see info

        Sources.ConstantMMF theta_PM(theta=PM.H_cB*l_PM)
        "Permanent magnet's magnetomotive force"
          annotation (extent=[-50,54; -30,74], rotation=270);
        MagneticGround magGround annotation (extent=[-50,28; -30,48]);
        Modelica_Magnetic.FluxTube.FixedShape.HollowCylinderRadialFlux R_mPM(
          my_rConst=PM.my_r,
          b=t,
          r_i=r_core + l_air,
          r_o=r_core + l_air + l_PM,
        nonLinearPermeability=false) "Reluctance of permanent magnet"
                                annotation (extent=[-32,70; -12,90]);
        Material.HardMagnetic.PermanentMagnetBehaviour PM(redeclare record
          material =   Modelica_Magnetic.Material.HardMagnetic.NdFeB, T_opCelsius=coil.T_opCelsius)
        "Permanent magnet material; coercitivity and relative permeabiliity used in theta_PM and R_mPM"
          annotation (extent=[-68,70; -48,90]);
        annotation (
        Images(Parameters(source="Images/Magnetic/Examples/ElectrodynamicActuator/MovingCoilActuator_dimensions.png")),
        Diagram,
          Documentation(info="<html>
<p>
Please refer to the <b>Parameters</b> section for a schematic drawing of this axisymmetric moving coil actuator. The half-section below shows the field lines of the permanent magnetic field (without armature current) obtained with a finite element analysis (FEA). The overlaid network of magnetic flux tube elements with nearly homogenous flux is rather simple in this example. Leakage fields are accounted for with a leakage coefficient and an appropriate leakage reluctance R_mLeak. Despite its simplicity, the model is well-suited for initial rough design of such moving coil actuators prior to detailed magnetic design, e.g. with FEA.
<br>
<br>
<dl>
<dd>
<p><IMG SRC=\"../Images/Magnetic/Examples/ElectrodynamicActuator/MagneticCircuitModel_fluxTubePartitioning.png\" ALT=\"Field lines and assigned flux tubes of the electrodynamic actuator\"></p>
</dd>
</dl>
<br>
During model-based actuator design, the radii of the flux tube elements (and hence their cross-sectional areas and flux densities) should be assigned with parametric equations so that common design rules are met (e.g. allowed flux density in ferromagnetic stator, operating point of permanent magnet on its demagnetization curve). For simplicity, those equations are omitted in the example. Instead, the found radii are assigned to the flux tube elements directly.
<br>
<br>
Only the permanent magnetic flux is considered in the lumped magnetic network. The influence of the coil-imposed magnetomotive force on the flux distribution in the magnetic subsystem is neglected. Also, the dependence of the coil inductance L_coil on the armature position x and thereof resulting reluctance forces are neglected. Consideration of these effects is possible in principle, but requires a higher modelling effort for the lumped magnetic network.
<br>
<br>
The coil inductance L_coil is calculated on base of an estimate for the total magnetic reluctance R_mTot that is \"seen\" by the coil. It is worth to note that the flux path for the coil-imposed magnetic flux is different from the flux path of the permanent magnetic flux. Only the portions of permanent magnet and air gap \"on the right side\" of the coil are part of the coil's main flux path. Although the formula for estimation of the total reluctance R_mTot is rather simple, comparison with FEA showed that it is well-suited for initial rough design and system simulation of the actuator. The relative difference to the inductance obtained with more accurate FEA is -12% for the armature in mid-position.
<br>
<br>
For useful operating currents, the relative differences of the air gap flux density and the resulting <i>Lorentz</i> force to the values obtained with FEA are within -5% to -8%. This accuracy is sufficient for initial actuator design and system simulation, too.
<br>
<br>
Have a look at <a href=\"FixedArmature\">FixedArmature</a> for an analysis of the actuator model under steady-state conditions. <a href=\"ArmatureStroke\">ArmatureStroke</a> is an example for the use of this model in dynamic simulation at system level.
<br>
<br>
</p>
</html>"),Icon(                 Rectangle(extent=[-90,100; 90,-100],   style(
              color=45,
              rgbcolor={255,128,0},
              fillColor=7,
              rgbfillColor={255,255,255})),
            Rectangle(extent=[-90,100; -50,-100], style(
                pattern=0,
                fillColor=45,
                rgbfillColor={255,128,0})),
            Rectangle(extent=[90,80; -90,100], style(
                pattern=0,
                fillColor=45,
                rgbfillColor={255,128,0},
                fillPattern=1)),
            Rectangle(extent=[90,-80; -90,-100], style(
                pattern=0,
                fillColor=45,
                rgbfillColor={255,128,0},
                fillPattern=1)),
            Rectangle(extent=[70,34; -90,-34], style(
                pattern=0,
                fillColor=45,
                rgbfillColor={255,128,0},
                fillPattern=1)),
            Rectangle(extent=[90,52; -16,64], style(
                pattern=0,
                fillColor=47,
                rgbfillColor={255,170,85})),
            Rectangle(extent=[90,-64; -12,-52], style(
                pattern=0,
                fillColor=47,
                rgbfillColor={255,170,85}))));
        Modelica_Magnetic.FluxTube.FixedShape.HollowCylinderRadialFlux R_mAir(
          b=t,
          my_rConst=1,
          r_i=r_core,
          r_o=r_core + l_air,
        nonLinearPermeability=false)
        "Reluctance of radial air gap between core and permanent magnet"
                         annotation (extent=[78,54; 98,74], rotation=270);
        Modelica_Magnetic.FluxTube.Leakage.LeakageWithCoefficient R_mLeak(
          c_leak=0.2,
          R_m(start=1e6),
          R_mUsefulTot=R_mAir.R_m + R_mFeCore.R_m + R_mFeBot.R_m + R_mFeOut.R_m)
        "Simple estimate for leakage reluctance"
                          annotation (extent=[-12,54; 8,74],   rotation=270);

        Modelica_Magnetic.FluxTube.FixedShape.HollowCylinderAxialFlux R_mFeCore(
          redeclare record Material =
            Modelica_Magnetic.Material.SoftMagnetic.Steel.Steel_9SMnPb28,
          l=t_add + t/2,
          r_o=r_core) "Reluctance of ferromagnetic stator core"
                   annotation (extent=[8,70; 28,90]);

        Modelica_Magnetic.FluxTube.FixedShape.HollowCylinderAxialFlux R_mFeOut(
          redeclare record Material =
            Modelica_Magnetic.Material.SoftMagnetic.Steel.Steel_9SMnPb28,
          l=t_add + t/2,
          r_i=r_core + l_air + l_PM,
          r_o=R_mFeOut.r_i + t_out)
        "Reluctance of outer section of ferromagnetic stator"
                       annotation (extent=[64,70; 84,90]);

        Modelica_Magnetic.Sources.CoilDesign coil(
          b_w=12e-3,
          U=6,
          J_desired=20e6,
          T_opCelsius=20,
          d_wireChosen=0.315e-3,
          w_chosen=140,
          h_w=l_air - 2*l_clear - l_car,
          l_avg=2*pi*(r_core + l_clear + l_car + coil.h_w/2))
        "Calculation of coil parameters (wire diameter, number of turns et al.) and recalculation with optionally chosen wire diameter"
                   annotation (extent=[48,20; 68,40]);
        Modelica.Electrical.Analog.Sensors.CurrentSensor iSensor "coil current"
          annotation (extent=[-90,0; -70,-20],   rotation=90);
        Modelica.Electrical.Analog.Sources.SignalVoltage u_i
        "induced voltage due to armature motion"
          annotation (extent=[-90,-50; -70,-30], rotation=270);
        Utilities.VariableGain c_force
        "Motor constant (gain variable due to variability of B_Air)"
          annotation (extent=[-50,-20; -30,0]);
        Utilities.VariableGain c_velocity
        "Motor constant (gain variable due to variability of B_Air)"
          annotation (extent=[-30,-50; -50,-30]);
        Modelica.Mechanics.Translational.Force MotorForce
          annotation (extent=[2,-20; 22,0]);
        Modelica.Mechanics.Translational.Sensors.SpeedSensor vSensor
        "armature velocity"
          annotation (extent=[30,-36; 50,-16],
                                             rotation=270);
        Utilities.TranslatoryArmature armature(
          m=m_a,
          x_max=x_max,
          x_min=x_min)
        "Inertia of armature and stoppers at end of stroke range"
                       annotation (extent=[60,-20; 80,0]);
        Modelica.Electrical.Analog.Basic.VariableInductor L_coil
          annotation (extent=[-90,10; -70,30], rotation=270);
        Modelica.Electrical.Analog.Basic.Resistor R_coil(R=coil.R_actual)
          annotation (extent=[-90,40; -70,60], rotation=270);
        Modelica_Magnetic.FluxTube.FixedShape.HollowCylinderRadialFlux R_mFeBot(
          redeclare record Material =
            Modelica_Magnetic.Material.SoftMagnetic.Steel.Steel_9SMnPb28,
          b=t_bot,
          r_i=r_core/2,
          r_o=R_mFeOut.r_i + t_out/2) annotation (extent=[36,70; 56,90]);

      equation
        c_force.k = R_mAir.B * coil.w_chosen*coil.l_avg;  // Motor constant
        c_velocity.k = c_force.k;
        R_mTot = 1/(R_mPM.G_m/2 + R_mLeak.G_m) + R_mFeCore.R_m + R_mFeBot.R_m + R_mFeOut.R_m + 2*R_mAir.R_m;   // see info
        L_coil.L = coil.w_chosen^2 / R_mTot;
        x = flange.s;

        connect(magGround.p, theta_PM.n) annotation (points=[-40,48; -40,54],  style(color=45, rgbcolor={
                255,127,0}));
        connect(theta_PM.p, R_mPM.p) annotation (points=[-40,74; -40,80; -32,80],
            style(color=45, rgbcolor={255,127,0}));
        connect(R_mLeak.p, R_mPM.n) annotation (points=[-2,74; -2,80; -12,80],
            style(color=45, rgbcolor={255,127,0}));
        connect(R_mAir.p, R_mFeOut.n)
                                    annotation (points=[88,74; 88,80; 84,80],
            style(color=45, rgbcolor={255,127,0}));
        connect(iSensor.n,u_i. p)   annotation (points=[-80,-20; -80,-30], style(
              color=3, rgbcolor={0,0,255}));
        connect(c_velocity.y,u_i. v) annotation (points=[-51,-40; -73,-40],
            style(
            color=74,
            rgbcolor={0,0,127},
            fillColor=0,
            rgbfillColor={0,0,0},
            fillPattern=1));
        connect(iSensor.i,c_force. u) annotation (points=[-70,-10; -52,-10],
            style(
            color=74,
            rgbcolor={0,0,127},
            fillColor=0,
            rgbfillColor={0,0,0},
            fillPattern=1));
        connect(u_i.n, n) annotation (points=[-80,-50; -80,-60; -100,-60],
            style(color=3, rgbcolor={0,0,255}));
        connect(MotorForce.flange_b,vSensor. flange_a) annotation (points=[22,-10;
              40,-10; 40,-16],
                     style(
            color=58,
            rgbcolor={0,127,0},
            fillPattern=1));
        connect(armature.flange_b, flange) annotation (points=[80,-10; 90,-10;
              90,0; 100,0],
            style(
            color=58,
            rgbcolor={0,127,0},
            pattern=0,
            fillColor=47,
            rgbfillColor={255,170,85},
            fillPattern=1));
        connect(vSensor.v, c_velocity.u) annotation (points=[40,-37; 40,-40; -28,-40],
                        style(
            color=74,
            rgbcolor={0,0,127},
            fillColor=0,
            rgbfillColor={0,0,0},
            fillPattern=1));
        connect(MotorForce.f, c_force.y) annotation (points=[0,-10; -29,-10],
            style(color=74, rgbcolor={0,0,127}));
        connect(L_coil.n, iSensor.p) annotation (points=[-80,10; -80,0], style(
              color=3, rgbcolor={0,0,255}));
        connect(R_coil.p, p)
          annotation (points=[-80,60; -100,60], style(color=3, rgbcolor={0,0,255}));
        connect(R_coil.n, L_coil.p)
          annotation (points=[-80,40; -80,30], style(color=3, rgbcolor={0,0,255}));
        connect(armature.flange_a, MotorForce.flange_b)
          annotation (points=[60,-10; 22,-10],
                                           style(color=58, rgbcolor={0,127,0}));
        connect(R_mLeak.n, theta_PM.n) annotation (points=[-2,54; -2,48; -40,48;
            -40,54],
            style(color=45, rgbcolor={255,127,0}));
        connect(R_mLeak.n, R_mAir.n) annotation (points=[-2,54; -2,48; 88,48;
            88,54],
            style(color=45, rgbcolor={255,127,0}));
        connect(R_mFeOut.p, R_mFeBot.n) annotation (points=[64,80; 56,80],
            style(color=45, rgbcolor={255,127,0}));
        connect(R_mFeBot.p, R_mFeCore.n) annotation (points=[36,80; 28,80],
            style(color=45, rgbcolor={255,127,0}));
        connect(R_mFeCore.p, R_mPM.n) annotation (points=[8,80; -12,80], style(
              color=45, rgbcolor={255,127,0}));
      end MagneticCircuitModel;

      model FixedArmature
      "MagneticCircuitModel with fixed armature; suited for analysis of actuator model during magnetic design"

        extends Modelica.Icons.Example;

        Modelica.Electrical.Analog.Sources.StepVoltage u_step(
            startTime=0, V=actuator.coil.U)
        "Supply voltage set to value defined in actuator.coil"
          annotation (extent=[-60,0; -40,20],  rotation=270);
        annotation (Diagram,                                       experiment(StopTime=
              0.1),
          Documentation(info="<html>
<p>
Have a look at <a href=\"Modelica_Magnetic.Examples.ElectrodynamicActuator\">ElectrodynamicActuator</a> for general comments and at <a href=\"Modelica_Magnetic.Examples.ElectrodynamicActuator.MagneticCircuitModel\">MagneticCircuitModel</a> for a detailed explanation of this actuator model. <br>
<br>
With the armature position fixed, all variables can easily be analyzed during rough design of the actuator. <b>Simulate for 0.1 s</b> and analyze, e.g. plot vs. time variables of interest such as flux densities:
<pre>
    actuator.R_mFeCore.B        // flux density in stator core
    actuator.R_mFeOut.B         // flux density in outer stator section
    actuator.R_mAir.B           // air gap flux density relevant for Lorentz force
    actuator.R_mPM.B            // flux density in permanent magnet.

</pre>
Compare the coil parameters calculated for desired operating conditions with the ones obtained for a chosen, available wire diameter:
<pre>
    actuator.coil.d_wireCalculated   vs.   actuator.coil.d_wireChosen           // wire diameter
    actuator.coil.w_calculated       vs.   actuator.coil.w_chosen               // number of turns
    actuator.coil.c_condFillChosen   vs.   actuator.coil.c_condFillActual       // conductor filling factor
    actuator.coil.J_desired          vs.   actuator.coil.J_actual               // current density.

</pre>
Watch the current rise due to the inductance of the armature coil by plotting vs. time:
<pre>
    actuator.p.i                 // input current to actuator.

</pre>
</p>
</html>"));
        Modelica.Electrical.Analog.Basic.Ground ground
          annotation (extent=[-60,-30; -40,-10]);
        MagneticCircuitModel actuator
          annotation (extent=[-20,0; 0,20]);
        Modelica.Mechanics.Translational.Fixed fixedPos(s0=(actuator.x_min + actuator.x_max)/2)
        "Fixed armature position"
          annotation (extent=[10,0; 30,20]);
      equation
        connect(ground.p, u_step.n)   annotation (points=[-50,-10; -50,0],
            style(color=3, rgbcolor={0,0,255}));
        connect(actuator.p, u_step.p)           annotation (points=[-20,16; -36,
              16; -36,20; -50,20],style(color=3, rgbcolor={0,0,255}));
        connect(actuator.n, u_step.n)           annotation (points=[-20,4; -36,
              4; -36,0; -50,0],      style(color=3, rgbcolor={0,0,255}));
        connect(fixedPos.flange_b, actuator.flange)
          annotation (points=[20,10; 0,10],
                                          style(color=58, rgbcolor={0,127,0}));
      end FixedArmature;

      model ArmatureStroke
      "Armature stroke of MagneticCircuitModel with load mass"

        extends Modelica.Icons.Example;

        Modelica.Electrical.Analog.Sources.StepVoltage u_step(V=actuator.coil.U,
            startTime=0)
          annotation (extent=[-60,0; -40,20],  rotation=270);
        annotation (Diagram,                                       experiment(StopTime=
                0.05),
          Documentation(info="<html>
<p>
Have a look at <a href=\"Modelica_Magnetic.Examples.ElectrodynamicActuator\">ElectrodynamicActuator</a> for general comments and at <a href=\"Modelica_Magnetic.Examples.ElectrodynamicActuator.MagneticCircuitModel\">MagneticCircuitModel</a> for a detailed explanation of this actuator model.
<br>
<br>
A voltage step is applied to the detailed model of the moving coil actuator. The actuator's armature and a therewith connected load mass perform a stroke between the two stoppers included in actuator.armature. <b>Simulate for 0.05 s</b> and plot vs. time:
<pre>
    actuator.p.i              // input current to actuator
    actuator.MotorForce.f     // <i>Lorentz</i> force, proportional to current
    m_load.v                  // armature and load mass velocity
    actuator.x                // armature and load mass position
</pre>
The initial current rise is due to the inductance of the actuator coil. After acceleration of armature and load, the actuator approaches to an equilibrium between the applied voltage and the motion-induced back-emf (idle operating conditions). Bouncing occurs when the armature arrives at the stopper actuator.armature.stopper_xMax. The bouncing is rather intense due to the absence of any kind of external friction in this simple example (apart from the nonlinear damping in the stopper elements). After decay of this bouncing, the actuator operates under conditions valid for a blocked armature.
</p>
</html>"));
        Modelica.Electrical.Analog.Basic.Ground ground
          annotation (extent=[-60,-30; -40,-10]);
        MagneticCircuitModel actuator
          annotation (extent=[-20,0; 0,20]);
        Modelica.Mechanics.Translational.SlidingMass m_load(m=0.05)
        "Load to be moved in addition to the armature mass"
                                 annotation (extent=[20,0; 40,20]);
      equation
        connect(ground.p, u_step.n)   annotation (points=[-50,-10; -50,0],
            style(color=3, rgbcolor={0,0,255}));
        connect(actuator.p, u_step.p)           annotation (points=[-20,16; -36,
              16; -36,20; -50,20],style(color=3, rgbcolor={0,0,255}));
        connect(actuator.n, u_step.n)           annotation (points=[-20,4; -36,
              4; -36,0; -50,0],      style(color=3, rgbcolor={0,0,255}));
        connect(m_load.flange_a, actuator.flange)
          annotation (points=[20,10; 0,10],
                                          style(color=58, rgbcolor={0,127,0}));
      end ArmatureStroke;

    end ElectrodynamicActuator;

    annotation (Documentation(info="<html>
</html>"));
    package ElectromagneticActuator
    "Two models of a reluctance actuator of different modelling depth and their comparison and usage"

      model SimpleSolenoidModel
      "Simple network model of a lifting magnet with planar armature end face"
        extends Modelica_Magnetic.Interfaces.ElectromechanicalActuator;

      //armature
        parameter SI.Radius r_arm = 5e-3 "Armature radius = pole radius";
        parameter SI.Length l_arm = 26e-3 "Armature length";

      //yoke
        parameter SI.Radius r_yokeOut = 15e-3 "Outer yoke radius";
        parameter SI.Radius r_yokeIn = 13.5e-3 "Inner yoke radius";
        parameter SI.Length l_yoke = 35e-3 "Axial yoke length";
        parameter SI.Thickness t_yokeBot = 3.5e-3
        "Axial thickness of yoke bottom";

      //pole
        parameter SI.Length l_pole = 6.5e-3 "Axial length of pole";
        parameter SI.Thickness t_poleBot = 3.5e-3
        "Axial thickness of bottom at pole side";

        parameter SI.Thickness t_airPar = 0.65e-3
        "Radial thickness of parasitic air gap due to slide guiding";

        parameter SI.Position x_max = 5e-3
        "Stopper at maximum armature position";
        parameter SI.Position x_min = 0.25e-3
        "Stopper at minimum armature position";

        SI.Position x(start=x_max, stateSelect=StateSelect.prefer)
        "Armature position, alias for flange position (identical with length of working air gap)";

    protected
        parameter SI.Density rho_steel = 7853
        "Density for calculation of armature mass from geometry";

        annotation (Images(Parameters(source="Images/Magnetic/Examples/ElectromagneticActuator/Solenoid_dimensions.png")),
          Coordsys(
            extent=[-100, -100; 100, 100],
            grid=[2, 2],
            component=[20, 20]),
          Window(
            x=0.16,
            y=0.15,
            width=0.6,
            height=0.6),
          Icon(
            Rectangle(extent=[90,-30; -4,30], style(
                color=45,
                rgbcolor={255,128,0},
                fillColor=45,
                rgbfillColor={255,128,0})),
            Rectangle(extent=[-40,-30; -90,30], style(
                color=45,
                rgbcolor={255,128,0},
                fillColor=45,
                rgbfillColor={255,128,0})),
            Rectangle(extent=[-80,-100; -90,100], style(
                color=45,
                rgbcolor={255,128,0},
                fillColor=45,
                rgbfillColor={255,128,0})),
            Rectangle(extent=[90,90; -90,100], style(
                color=45,
                rgbcolor={255,128,0},
                fillColor=45,
                rgbfillColor={255,128,0})),
            Rectangle(extent=[90,-100; -90,-90], style(
                color=45,
                rgbcolor={255,128,0},
                fillColor=45,
                rgbfillColor={255,128,0})),
            Rectangle(extent=[90,40; 80,100], style(
                color=45,
                rgbcolor={255,128,0},
                fillColor=45,
                rgbfillColor={255,128,0})),
            Rectangle(extent=[90,-100; 80,-40], style(
                color=45,
                rgbcolor={255,128,0},
                fillColor=45,
                rgbfillColor={255,128,0})),
            Rectangle(extent=[-70,80; 70,40], style(
                color=52,
                rgbcolor={255,213,170},
                fillColor=52,
                rgbfillColor={255,213,170},
                fillPattern=1)),
            Rectangle(extent=[-70,-40; 70,-80], style(
                color=52,
                rgbcolor={255,213,170},
                fillColor=52,
                rgbfillColor={255,213,170},
                fillPattern=1))),
          Diagram,
        Documentation(info="<HTML>
<p>
Please refer to the <b>Parameters</b> section for a schematic drawing of this axisymmetric lifting magnet.
In the half-section below, the flux tube elements of the actuator's magnetic circuit are superimposed on a field plot obtained with FEA. The magnetomotive force imposed by the coil is modelled as one lumped element. As a result, the radial leakage flux between armature and yoke that occurs especially at large working air gaps can not be considered properly. This leads to a a higher total reluctance and lower inductance respectively compared to FEA for large working air gaps (i.e. armature close to x_max). Please have a look at the comments associated with the individual model components for a short explanation of their purpose in the model.<p>
<br>
</p>
<p>
<dl>
<dd>
<p><IMG SRC=\"../Images/Magnetic/Examples/ElectromagneticActuator/SimpleSolenoidModel_fluxTubePartitioning.png\" ALT=\"Field lines and assigned flux tubes of the simple solenoid model\"></p>
</dd>
</dl>
</p>
<p>
The coupling coefficient c_coupl in the coil is set to 1 in this example, since leakage flux is accounted for explicitly with the flux tube element G_mLeakWork. Although this leakage model is rather simple, it describes the reluctance force due to the leakage field sufficiently, especially at large air gaps. With decreasing air gap length, the influence of the leakage flux on the actuator's net reluctance force decreases due to the increasing influence of the main working air gap G_mAirWork.
</p>
<p>
During model-based actuator design, the radii and lengths of the flux tube elements (and hence their cross-sectional areas and flux densities) should be assigned with parametric equations so that common design rules are met (e.g. allowed flux density in ferromagnetic parts, allowed current density and required cross-sectional area of winding). For simplicity, those equations are omitted in the example. Instead, the found values are assigned to the model elements directly.
</p>
</HTML>"));

    public
        MagneticGround magGround      annotation (extent=[46,4; 66,24],  rotation=0);
        Sources.ElectroMagneticConverter coil(           c_coupl=1, w=957)
        "Electro-magnetic converter"
          annotation (extent=[-14,10; 6,30],  rotation=90);
        Modelica.Electrical.Analog.Basic.Resistor R_coil(R=10)
        "Coil resistance"
          annotation (extent=[-60,-30; -40,-10],
                                               rotation=0);
        Modelica_Magnetic.FluxTube.FixedShape.HollowCylinderAxialFlux
        G_mFeYokeSide(
          redeclare record Material =
            Modelica_Magnetic.Material.SoftMagnetic.Steel.Steel_9SMnPb28,
          l=l_yoke - (t_poleBot + t_yokeBot)/2,
          r_i=r_yokeIn,
          r_o=r_yokeOut)
        "Permeance of of hollow cylindric section of ferromagnetic yoke"
          annotation (extent=[-10,70; 10,90]);

        Modelica_Magnetic.FluxTube.FixedShape.HollowCylinderAxialFlux G_mFeArm(
          r_i=0,
          redeclare record Material =
            Modelica_Magnetic.Material.SoftMagnetic.Steel.Steel_9SMnPb28,
          l=l_yoke - (t_yokeBot + t_poleBot)/2 - l_pole - (x_max + x_min)/2,
          r_o=r_arm) "Permeance of ferfomagnetic armature"
          annotation (extent=[8,20; 28,40], rotation=180);

        FluxTube.Force.HollowCylinderAxialFlux G_mAirWork(r_o=r_arm)
        "Permeance of working air gap (between armature and pole end faces)"
          annotation (extent=[-42,20; -22,40],
                                             rotation=180);
        Modelica_Magnetic.FluxTube.FixedShape.HollowCylinderRadialFlux
        G_mFeYokeBot(
          redeclare record Material =
            Modelica_Magnetic.Material.SoftMagnetic.Steel.Steel_9SMnPb28,
          b=t_yokeBot,
          r_i=r_arm + t_airPar,
          r_o=r_yokeIn) "Permeance of bottom side of ferromagnetic yoke"
          annotation (extent=[60,54; 80,74],   rotation=90);

        Modelica_Magnetic.FluxTube.FixedShape.HollowCylinderRadialFlux
        G_mAirPar(
          my_rConst=1,
          b=t_yokeBot,
          r_i=r_arm,
          r_o=r_arm + t_airPar,
        nonLinearPermeability=false)
        "Permeance of parasitic radial air gap due to slide guiding"
          annotation (extent=[60,30; 80,50],   rotation=90);
        Modelica_Magnetic.FluxTube.FixedShape.HollowCylinderRadialFlux
        G_mFePoleBot(
          redeclare record Material =
            Modelica_Magnetic.Material.SoftMagnetic.Steel.Steel_9SMnPb28,
          b=t_poleBot,
          r_i=r_arm,
          r_o=r_yokeIn) "Permeance of bottom side of pole"
          annotation (extent=[-82,54; -62,74],
                                             rotation=90);

        Modelica_Magnetic.FluxTube.FixedShape.HollowCylinderAxialFlux G_mFePole(
          redeclare record Material =
            Modelica_Magnetic.Material.SoftMagnetic.Steel.Steel_9SMnPb28,
          l=l_pole,
          r_o=r_arm) "Permeance of ferromagnetic pole"
          annotation (extent=[-70,40; -50,20],
                                             rotation=180);

        Modelica_Magnetic.Utilities.TranslatoryArmature armature(
          x_max=x_max,
          x_min=x_min,
          m=rho_steel*l_arm*pi*r_arm^2)
        "Inertia of armature and stoppers at end of stroke range"
                                    annotation (extent=[64,-10; 84,10]);
        FluxTube.Leakage.QuarterCylinder G_mLeak1(t=2*pi*(r_arm + t_airPar/2))
        "Leakage pereance between inner edge of yoke bore and armature side face"
          annotation (extent=[46,30; 66,50],   rotation=90);
        FluxTube.Leakage.QuarterHollowCylinder G_mLeak2(            ratio=8, t=2*pi*
              r_arm)
        "Leakage permeance between inner side of yoke bottom and armature side (r_i = t_airPar)"
                           annotation (extent=[32,30; 52,50],   rotation=90);
        FluxTube.Force.LeakageAroundPoles G_mLeakWork(                           t=2*pi*
              r_arm, r_leak=0.003)
        "Permeance of leakage air gap around working air gap (between armature and pole side faces)"
          annotation (extent=[-42,54; -22,34],
                                             rotation=180);
      equation
        x = flange.s;
        connect(R_coil.p, p) annotation (points=[-60,-20; -90,-20; -90,60; -100,60],
                                                                           style(
              color=3, rgbcolor={0,0,255}));
        connect(R_coil.n, coil.p_el) annotation (points=[-40,-20; -10,-20; -10,10],
                                                                                style(
              color=3, rgbcolor={0,0,255}));
        connect(coil.n_el, n) annotation (points=[2,10; 2,-60; -100,-60],    style(
              color=3, rgbcolor={0,0,255}));
        connect(armature.flange_b, flange)     annotation (points=[84,0; 100,0], style(
              color=58, rgbcolor={0,127,0}));
        connect(armature.flange_a, G_mAirWork.flange) annotation (points=[64,0; -32,0;
              -32,30],      style(color=58, rgbcolor={0,127,0}));
        connect(G_mAirWork.flange, G_mLeakWork.flange)
                                                   annotation (points=[-32,30; -32,44],
                      style(color=58, rgbcolor={0,127,0}));
        connect(G_mFeYokeBot.n, G_mFeYokeSide.n) annotation (points=[70,74; 70,80; 10,
              80], style(color=45, rgbcolor={255,127,0}));
        connect(G_mFePoleBot.n, G_mFeYokeSide.p) annotation (points=[-72,74; -72,80;
              -10,80], style(color=45, rgbcolor={255,127,0}));
        connect(G_mLeakWork.n, G_mAirWork.n) annotation (points=[-42,44; -42,30],
            style(color=45, rgbcolor={255,127,0}));
        connect(G_mAirWork.p, G_mLeakWork.p) annotation (points=[-22,30; -22,44],
            style(color=45, rgbcolor={255,127,0}));
        connect(G_mAirWork.n, G_mFePole.p) annotation (points=[-42,30; -50,30], style(
              color=45, rgbcolor={255,127,0}));
        connect(G_mFePole.n, G_mFePoleBot.p) annotation (points=[-70,30; -72,30; -72,
              54], style(color=45, rgbcolor={255,127,0}));
        connect(coil.n_mag, G_mFeArm.n)
          annotation (points=[2,30; 8,30], style(color=45, rgbcolor={255,127,0}));
        connect(G_mFeArm.p, G_mLeak2.p)
          annotation (points=[28,30; 42,30], style(color=45, rgbcolor={255,127,0}));
        connect(G_mLeak2.p, G_mLeak1.p)
          annotation (points=[42,30; 56,30], style(color=45, rgbcolor={255,127,0}));
        connect(G_mLeak1.p, G_mAirPar.p)
          annotation (points=[56,30; 70,30], style(color=45, rgbcolor={255,127,0}));
        connect(magGround.p, G_mLeak1.p) annotation (points=[56,24; 56,30],
            style(color=45, rgbcolor={255,127,0}));
        connect(G_mAirPar.n, G_mLeak1.n)
          annotation (points=[70,50; 56,50], style(color=45, rgbcolor={255,127,0}));
        connect(G_mLeak1.n, G_mLeak2.n)
          annotation (points=[56,50; 42,50], style(color=45, rgbcolor={255,127,0}));
        connect(G_mAirPar.n, G_mFeYokeBot.p)
          annotation (points=[70,50; 70,54], style(color=45, rgbcolor={255,127,0}));
        connect(coil.p_mag, G_mAirWork.p) annotation (points=[-10,30; -22,30], style(
              color=45, rgbcolor={255,127,0}));
      end SimpleSolenoidModel;

      model AdvancedSolenoidModel
      "Advanced network model of a lifting magnet with planar armature end face, split magnetomotive force"
        extends Modelica_Magnetic.Interfaces.ElectromechanicalActuator;

      //armature
        parameter SI.Radius r_arm = 5e-3 "Armature radius = pole radius";
        parameter SI.Length l_arm = 26e-3 "Armature length";

      //yoke
        parameter SI.Radius r_yokeOut = 15e-3 "Outer yoke radius";
        parameter SI.Radius r_yokeIn = 13.5e-3 "Inner yoke radius";
        parameter SI.Length l_yoke = 35e-3 "Axial yoke length";
        parameter SI.Thickness t_yokeBot = 3.5e-3
        "Axial thickness of yoke bottom";

      //pole
        parameter SI.Length l_pole = 6.5e-3 "Axial length of pole";
        parameter SI.Thickness t_poleBot = 3.5e-3
        "Axial thickness of bottom at pole side";

        parameter SI.Thickness t_airPar = 0.65e-3
        "Radial thickness of parasitic air gap due to slide guiding";

        parameter SI.Position x_max = 5e-3
        "Stopper at maximum armature position";
        parameter SI.Position x_min = 0.25e-3
        "Stopper at minimum armature position";

        SI.Position x(start=x_max, stateSelect=StateSelect.prefer)
        "Armature position";

        parameter Real w = 957 "Number of turns";

        SI.MagneticFlux Psi_tot "Total flux linkage for information only";
        SI.Inductance L_statTot "Total static inductance for information only";

    protected
        parameter SI.Density rho_steel = 7853
        "Density for calculation of armature mass from geometry";

        annotation (Images(Parameters(source="Images/Magnetic/Examples/ElectromagneticActuator/Solenoid_dimensions.png")),
          Coordsys(
            extent=[-100, -100; 100, 100],
            grid=[2, 2],
            component=[20, 20]),
          Window(
            x=0.16,
            y=0.15,
            width=0.6,
            height=0.6),
          Icon(
            Rectangle(extent=[90,-30; -4,30], style(
                color=45,
                rgbcolor={255,128,0},
                fillColor=45,
                rgbfillColor={255,128,0})),
            Rectangle(extent=[-40,-30; -90,30], style(
                color=45,
                rgbcolor={255,128,0},
                fillColor=45,
                rgbfillColor={255,128,0})),
            Rectangle(extent=[-80,-100; -90,100], style(
                color=45,
                rgbcolor={255,128,0},
                fillColor=45,
                rgbfillColor={255,128,0})),
            Rectangle(extent=[90,90; -90,100], style(
                color=45,
                rgbcolor={255,128,0},
                fillColor=45,
                rgbfillColor={255,128,0})),
            Rectangle(extent=[90,-100; -90,-90], style(
                color=45,
                rgbcolor={255,128,0},
                fillColor=45,
                rgbfillColor={255,128,0})),
            Rectangle(extent=[90,40; 80,100], style(
                color=45,
                rgbcolor={255,128,0},
                fillColor=45,
                rgbfillColor={255,128,0})),
            Rectangle(extent=[90,-100; 80,-40], style(
                color=45,
                rgbcolor={255,128,0},
                fillColor=45,
                rgbfillColor={255,128,0})),
            Rectangle(extent=[-70,80; 70,40], style(
                color=52,
                rgbcolor={255,213,170},
                fillColor=52,
                rgbfillColor={255,213,170},
                fillPattern=1)),
            Rectangle(extent=[-70,-40; 70,-80], style(
                color=52,
                rgbcolor={255,213,170},
                fillColor=52,
                rgbfillColor={255,213,170},
                fillPattern=1)),
            Line(points=[4,30; 4,32; 2,38; -4,48; -14,60; -22,72; -24,80; -24,
                  90], style(color=45, rgbcolor={255,128,0})),
            Line(points=[22,30; 22,32; 20,38; 14,48; 4,60; -4,72; -6,80; -6,90],
                style(color=45, rgbcolor={255,128,0})),
            Line(points=[40,30; 40,32; 38,38; 32,48; 22,60; 14,72; 12,80; 12,90],
                style(color=45, rgbcolor={255,128,0}))),
          Diagram,
        Documentation(info="<HTML>
<p>
Please have a look at <a href=\"SimpleSolenoidModel\">SimpleSolenoidModel</a> for a general description of this actuator. Unlike in that simple magnetic network model, the coil is split into two lumped elements here. This enables for more realistic modelling of the radial leakage flux between armature and yoke (leakage permeance G_mLeakRad). Especially for large air gaps, the influence of this leakage flux on the actuator's inductance and its electromagnetic force is rather strong. Please have a look at <a href=\"ComparisonQuasiStationary\">ComparisonQuasiStationary</a> for a comparison of both models with FEA-based results included as reference.
<br>
<br>
<dl>
<dd>
<p><IMG SRC=\"../Images/Magnetic/Examples/ElectromagneticActuator/AdvancedSolenoidModel_fluxTubePartitioning.png\" ALT=\"Field lines and assigned flux tubes of the advanced solenoid model\"></p>
</dd>
</dl>
<br>
<p>
The parasitic capacitances c_par1 and c_par2 accross both partial coils assure that the voltages across these coils are well-defined during simulation.
</p>
</HTML>"));

    public
        MagneticGround magGround      annotation (extent=[42,2; 62,22],  rotation=0);
        Sources.ElectroMagneticConverter coil1(c_coupl=1, w=w/2)
        "Electro-magnetic conversion in first half of coil"
          annotation (extent=[-56,10; -36,30],rotation=90);
        Modelica.Electrical.Analog.Basic.Resistor R_coil1(R=5)
        "Resistance of first half of coil"
          annotation (extent=[-84,-30; -64,-10],
                                               rotation=0);
        Modelica_Magnetic.FluxTube.FixedShape.HollowCylinderAxialFlux
        G_mFeYokeSide1(
          redeclare record Material =
            Modelica_Magnetic.Material.SoftMagnetic.Steel.Steel_9SMnPb28,
          l=l_yoke/2 - t_poleBot/2,
          r_i=r_yokeIn,
          r_o=r_yokeOut)
        "Permeance of of first half of yoke's hollow cylindric section"
          annotation (extent=[-50,70; -30,90]);

        Modelica_Magnetic.FluxTube.FixedShape.HollowCylinderAxialFlux G_mFeArm(
          r_i=0,
          redeclare record Material =
            Modelica_Magnetic.Material.SoftMagnetic.Steel.Steel_9SMnPb28,
          l=l_yoke - (t_yokeBot + t_poleBot)/2 - l_pole - (x_max + x_min)/2,
          r_o=r_arm) "Permeance of ferfomagnetic armature"
          annotation (extent=[0,40; 20,20], rotation=180);

        FluxTube.Force.HollowCylinderAxialFlux G_mAirWork(r_o=r_arm)
        "Permeance of working air gap (between armature and pole end faces)"
          annotation (extent=[-32,20; -12,40],
                                             rotation=180);
        Modelica_Magnetic.FluxTube.FixedShape.HollowCylinderRadialFlux
        G_mFeYokeBot(
          redeclare record Material =
            Modelica_Magnetic.Material.SoftMagnetic.Steel.Steel_9SMnPb28,
          b=t_yokeBot,
          r_i=r_arm + t_airPar,
          r_o=r_yokeIn) "Permeance of bottom side of ferromagnetic yoke"
          annotation (extent=[64,54; 84,74],   rotation=90);

        Modelica_Magnetic.FluxTube.FixedShape.HollowCylinderRadialFlux
        G_mAirPar(
          my_rConst=1,
          b=t_yokeBot,
          r_i=r_arm,
          r_o=r_arm + t_airPar,
          nonLinearPermeability=false)
        "Permeance of parasitic radial air gap due to slide guiding"
          annotation (extent=[64,30; 84,50],   rotation=90);
        Modelica_Magnetic.FluxTube.FixedShape.HollowCylinderRadialFlux
        G_mFePoleBot(
          redeclare record Material =
            Modelica_Magnetic.Material.SoftMagnetic.Steel.Steel_9SMnPb28,
          b=t_poleBot,
          r_i=r_arm,
          r_o=r_yokeIn) "Permeance of bottom side of pole"
          annotation (extent=[-88,46; -68,66],
                                             rotation=90);

        Modelica_Magnetic.FluxTube.FixedShape.HollowCylinderAxialFlux G_mFePole(
          redeclare record Material =
            Modelica_Magnetic.Material.SoftMagnetic.Steel.Steel_9SMnPb28,
          l=l_pole,
          r_o=r_arm) "Permeance of ferromagnetic pole"
          annotation (extent=[-78,40; -58,20],
                                             rotation=180);

        Modelica_Magnetic.Utilities.TranslatoryArmature armature(
          x_max=x_max,
          x_min=x_min,
          m=rho_steel*l_arm*pi*r_arm^2)
        "Inertia of armature and stoppers at end of stroke range"
                                    annotation (extent=[62,-10; 82,10]);
        FluxTube.Leakage.QuarterCylinder G_mLeak1(t=2*pi*(r_arm + t_airPar/2))
        "Leakage pereance between inner edge of yoke bore and armature side face"
          annotation (extent=[50,30; 70,50],   rotation=90);
        FluxTube.Leakage.QuarterHollowCylinder G_mLeak2(            ratio=8, t=2*pi*
              r_arm)
        "Leakage permeance between inner side of yoke bottom and armature side (r_i = t_airPar)"
                           annotation (extent=[36,30; 56,50],   rotation=90);
        Sources.ElectroMagneticConverter coil2(c_coupl=1, w=w/2)
        "Electro-magnetic conversion in first half of coil"
          annotation (extent=[20,10; 40,30],  rotation=90);
        Modelica.Electrical.Analog.Basic.Capacitor C_par1(C=1e-9)
        "Parasitic capacitance assigned to first half of coil"
          annotation (extent=[-60,-50; -40,-30], rotation=0);
        Modelica_Magnetic.FluxTube.FixedShape.HollowCylinderRadialFlux
        G_mLeakRad(
          my_rConst=1,
          r_i=r_arm,
          r_o=r_yokeIn,
          b=l_yoke/4,
        nonLinearPermeability=false)
        "Permeance of radial leakage flux tube between armature side and yoke side"
          annotation (extent=[-10,46; 10,66],  rotation=90);
        Modelica_Magnetic.FluxTube.FixedShape.HollowCylinderAxialFlux
        G_mFeYokeSide2(
          redeclare record Material =
            Modelica_Magnetic.Material.SoftMagnetic.Steel.Steel_9SMnPb28,
          l=l_yoke/2 - t_yokeBot/2,
          r_i=r_yokeIn,
          r_o=r_yokeOut)
        "Permeance of of second half of yoke's hollow cylindric section"
          annotation (extent=[20,70; 40,90]);

        Modelica.Electrical.Analog.Basic.Capacitor C_par2(C=1e-9)
        "Parasitic capacitance assigned to second half of coil"
          annotation (extent=[16,-50; 36,-30], rotation=0);
        Modelica.Electrical.Analog.Basic.Resistor R_par1(R=1e5)
        "Parasitic resistance assigned to first half of coil"
          annotation (extent=[-84,-50; -64,-30],
                                               rotation=0);
        Modelica.Electrical.Analog.Basic.Resistor R_par2(R=1e5)
        "Parasitic resistance assigned to second half of coil"
          annotation (extent=[-8,-50; 12,-30], rotation=0);
        Modelica.Electrical.Analog.Basic.Resistor R_coil2(R=5)
        "Resistance of second half of coil"
          annotation (extent=[-8,-30; 12,-10], rotation=0);
        FluxTube.Leakage.QuarterCylinder G_mLeak3(t=2*pi*(r_arm + t_airPar/2))
        "Leakage pereance between outer edge of yoke bore and armature side face"
          annotation (extent=[78,30; 98,50],   rotation=90);
        FluxTube.Force.LeakageAroundPoles G_mLeakWork(            r_leak=0.003, t=2*pi*
              r_arm)
        "Permeance of leakage air gap around working air gap (between armature and pole side faces)"
          annotation (extent=[-32,54; -12,34],
                                             rotation=180);
      equation
        x = flange.s;
        Psi_tot = coil1.Psi + coil2.Psi;
        L_statTot = coil1.L_stat + coil2.L_stat;
        connect(R_coil1.n, coil1.p_el)
                                     annotation (points=[-64,-20; -52,-20; -52,
              10],                                                              style(
              color=3, rgbcolor={0,0,255}));
        connect(armature.flange_b, flange)     annotation (points=[82,0; 100,0], style(
              color=58, rgbcolor={0,127,0}));
        connect(R_par1.n, C_par1.p)
                                 annotation (points=[-64,-40; -60,-40], style(
              color=3, rgbcolor={0,0,255}));
        connect(R_par1.p, R_coil1.p)  annotation (points=[-84,-40; -84,-20],
            style(color=3, rgbcolor={0,0,255}));
        connect(C_par2.p, R_par2.n)
                                 annotation (points=[16,-40; 12,-40],  style(
              color=3, rgbcolor={0,0,255}));
        connect(coil1.n_el,R_coil2. p) annotation (points=[-40,10; -40,-20; -8,
              -20],                                                       style(
              color=3, rgbcolor={0,0,255}));
        connect(R_coil2.n, coil2.p_el) annotation (points=[12,-20; 24,-20; 24,
              10],style(color=3, rgbcolor={0,0,255}));
        connect(R_par2.p, R_coil2.p)  annotation (points=[-8,-40; -8,-20],
            style(color=3, rgbcolor={0,0,255}));
        connect(R_coil1.p, p) annotation (points=[-84,-20; -92,-20; -92,60;
              -100,60],
            style(color=3, rgbcolor={0,0,255}));
        connect(coil1.p_mag, G_mFePole.p) annotation (points=[-52,30; -58,30],
            style(color=45, rgbcolor={255,127,0}));
        connect(G_mAirWork.n, coil1.n_mag) annotation (points=[-32,30; -40,30],
            style(color=45, rgbcolor={255,127,0}));
        connect(G_mAirWork.p, G_mLeakWork.p) annotation (points=[-12,30; -12,44],
            style(color=45, rgbcolor={255,127,0}));
        connect(G_mAirWork.n, G_mLeakWork.n) annotation (points=[-32,30; -32,44],
            style(color=45, rgbcolor={255,127,0}));
        connect(G_mLeakWork.flange, G_mAirWork.flange) annotation (points=[-22,44;
              -22,30],     style(color=58, rgbcolor={0,127,0}));
        connect(G_mAirWork.flange, armature.flange_a) annotation (points=[-22,30;
              -22,0; 62,0],       style(color=58, rgbcolor={0,127,0}));
        connect(n, C_par2.n) annotation (points=[-100,-60; 36,-60; 36,-40],
            style(color=3, rgbcolor={0,0,255}));
        connect(G_mFePole.n, G_mFePoleBot.p) annotation (points=[-78,30; -78,46],
            style(color=45, rgbcolor={255,127,0}));
        connect(G_mFePoleBot.n, G_mFeYokeSide1.p) annotation (points=[-78,66;
              -78,80; -50,80], style(color=45, rgbcolor={255,127,0}));
        connect(G_mFeYokeSide1.n, G_mFeYokeSide2.p) annotation (points=[-30,80;
              20,80], style(color=45, rgbcolor={255,127,0}));
        connect(G_mFeYokeSide2.n, G_mFeYokeBot.n) annotation (points=[40,80; 74,
              80; 74,74],
                   style(color=45, rgbcolor={255,127,0}));
        connect(G_mFeYokeBot.p, G_mAirPar.n) annotation (points=[74,54; 74,50],
            style(color=45, rgbcolor={255,127,0}));
        connect(G_mLeak3.n, G_mAirPar.n) annotation (points=[88,50; 74,50],
            style(color=45, rgbcolor={255,127,0}));
        connect(G_mAirPar.n, G_mLeak1.n) annotation (points=[74,50; 60,50],
            style(color=45, rgbcolor={255,127,0}));
        connect(G_mLeak1.n, G_mLeak2.n) annotation (points=[60,50; 46,50],
            style(color=45, rgbcolor={255,127,0}));
        connect(G_mLeak3.p, G_mAirPar.p) annotation (points=[88,30; 74,30],
            style(color=45, rgbcolor={255,127,0}));
        connect(G_mAirPar.p, G_mLeak1.p) annotation (points=[74,30; 60,30],
            style(color=45, rgbcolor={255,127,0}));
        connect(G_mLeak1.p, G_mLeak2.p) annotation (points=[60,30; 46,30],
            style(color=45, rgbcolor={255,127,0}));
        connect(coil2.n_mag, G_mLeak2.p) annotation (points=[36,30; 46,30],
            style(color=45, rgbcolor={255,127,0}));
        connect(coil2.p_mag, G_mFeArm.p) annotation (points=[24,30; 20,30],
            style(color=45, rgbcolor={255,127,0}));
        connect(G_mFeArm.n, G_mAirWork.p) annotation (points=[0,30; -12,30],
            style(color=45, rgbcolor={255,127,0}));
        connect(G_mFeArm.n, G_mLeakRad.p) annotation (points=[0,30; 0,46;
            -6.12303e-016,46],
                       style(color=45, rgbcolor={255,127,0}));
        connect(G_mLeakRad.n, G_mFeYokeSide1.n) annotation (points=[
            6.12303e-016,66; 6.12303e-016,80; -30,80],
                           style(color=45, rgbcolor={255,127,0}));
        connect(magGround.p, G_mLeak1.p) annotation (points=[52,22; 52,30; 60,
              30],
            style(color=45, rgbcolor={255,127,0}));
        connect(C_par2.n, coil2.n_el) annotation (points=[36,-40; 36,10], style(
              color=3, rgbcolor={0,0,255}));
        connect(C_par1.n, coil1.n_el) annotation (points=[-40,-40; -40,10],
            style(color=3, rgbcolor={0,0,255}));
      end AdvancedSolenoidModel;

      model ComparisonQuasiStationary
      "Slow forced armature motion of both solenoid models so that electromagnetic field and current are quasi-stationary"

        extends Modelica.Icons.Example;

        parameter SI.Voltage u_step = 12 "Applied voltage";

        Modelica.Electrical.Analog.Sources.StepVoltage u_source(V=u_step)
          annotation (extent=[-80,40; -60,60], rotation=270);
        annotation (Diagram,                                       experiment(StopTime=
                10, Tolerance=1e-007),
          Documentation(info="<html>
<p>
Have a look at <a href=\"ElectromagneticActuator\">ElectromagneticActuator</a> for general comments and at <a href=\"SimpleSolenoidModel\">SimpleSolenoidModel</a> and <a href=\"AdvancedSolenoidModel\">AdvancedSolenoidModel</a> for a detailed description of both magnetic network models.
</p>
<p>
Similar to static force-stroke measurements on real actuators, the armatures of both actuator models are forced to move slowly here. Hence, the dynamics of the electrical subsystems due to coil inductance and armature motion can be neglected and the static force-stroke characteristics are obtained. To illustrate the accuracy to be expected from the lumped magnetic network models, results obtained with stationary FEA are included as reference (position-dependent force, armature flux and actuator inductance). Note that these reference values are valid for the default supply voltage u_step=12V DC only!
</p>
<p>
Set the <b>tolerance</b> to <b>1e-7</b> and <b>simulate for 10 s</b>. Plot in one common window the electromagnetic force of the two magnetic network models and the FEA reference <b>vs. armature position x_set.y</b>:</p>
<pre>
    simpleMagnet.armature.flange_a.f      // electromagnetic force of simple magnetic network model
    advancedMagnet.armature.flange_a.f    // electromagnetic force of advaned magnetic network model
    comparisonWithFEA.y[1]                // electromagnetic force obtained with FEA as reference
</pre>
<p>
Electromagnetic or reluctance forces always act towards a decrease of air gap lengths. With the defined armature position coordinate x, the forces of the models are negative. </p>
<p>
The magnetic flux through the armature and the actuator's static inductance both illustrate the differences between the two magnetic network models. Similar to the forces, compare these quantities in one common plot window for each variable (plot vs. armature position x_set.y):</p>
<pre>
    simpleMagnet.G_mFeArm.Phi             // magnetic flux through armature of simple magnetic network model
    advancedMagnet.G_mFeArm.Phi           // magnetic flux through armature of advanced magnetic network model
    comparisonWithFEA.y[2]                // magnetic flux obtained with FEA as reference


    simpleMagnet.coil.L_stat              // static inductance of simple magnetic network model
    advancedMagnet.L_statTot              // series connection of both partial coils of advanced network model
    comparisonWithFEA.y[3]                // static inductance obtained with FEA as reference
</pre>
<p>
As mentioned in the description of both magnetic network models, one can tell the higher armature flux and inductance of the advanced solenoid model at large air gaps compared to that of the simple model. The effect of this difference on dynamic model behaviour can be analysed in <a href=\"ComparisonPullInStroke\">ComparisonPullInStroke</a>.
</p>
</html>"));
        Modelica.Electrical.Analog.Basic.Ground ground
          annotation (extent=[-80,10; -60,30]);
        Modelica.Mechanics.Translational.Position feed_x(f_crit=1000, exact=
              true)
                   annotation (extent=[-10,40; 10,60],rotation=180);
        Modelica.Blocks.Sources.Ramp x_set(
          duration=10,
          height=-(advancedMagnet.x_max - advancedMagnet.x_min),
          offset=advancedMagnet.x_max)
        "Prescribed armature position, slow enforced motion from x_max to x_min"
                               annotation (extent=[80,-10; 60,10]);
        Modelica.Mechanics.Translational.Position feed_x1(f_crit=1000, exact=
              false)
                   annotation (extent=[-10,-60; 10,-40],
                                                      rotation=180);
        Modelica.Electrical.Analog.Sources.StepVoltage u_source1(V=u_step)
          annotation (extent=[-80,-60; -60,-40],
                                               rotation=270);
        Modelica_Magnetic.Examples.ElectromagneticActuator.SimpleSolenoidModel
        simpleMagnet
          annotation (extent=[-40,-60; -20,-40]);
        Modelica_Magnetic.Examples.ElectromagneticActuator.AdvancedSolenoidModel
        advancedMagnet
          annotation (extent=[-40,40; -20,60]);
        Modelica.Blocks.Tables.CombiTable1Ds comparisonWithFEA(table=[0.00025,-85.8619,
              0.00014821,0.11954; 0.0005,-59.9662,0.00013931,0.11004; 0.00075,-41.0806,
              0.0001277,0.098942; 0.001,-28.88,0.00011587,0.088425; 0.00125,-21.4113,
              0.00010643,0.08015; 0.0015,-16.8003,9.9406e-005,0.073992; 0.00175,
              -13.6942,9.3416e-005,0.068792; 0.002,-11.1188,8.8564e-005,
              0.064492; 0.00225,-9.6603,8.4505e-005,0.060917; 0.0025,-8.4835,
              8.1215e-005,0.058017; 0.00275,-7.4658,7.7881e-005,0.055125; 0.003,
              -6.5591,7.5197e-005,0.052733; 0.00325,-5.9706,7.2447e-005,0.05035;
              0.0035,-5.5013,7.0342e-005,0.048525; 0.00375,-5.0469,6.8527e-005,
              0.046867; 0.004,-4.6573,6.6526e-005,0.045158; 0.00425,-4.2977,
              6.4425e-005,0.043442; 0.0045,-4.0912,6.2747e-005,0.04205; 0.00475,
              -3.7456,6.1231e-005,0.040733; 0.005,-3.5869,5.9691e-005,0.039467])
        "Valid for u_source=12V only; column 1: position, col.2: force, col.3: armature flux, col.4: inductance"
          annotation (extent=[60,60; 80,80]);
        Modelica.Electrical.Analog.Basic.Ground ground1
          annotation (extent=[-80,-90; -60,-70]);
      equation
        connect(ground.p, u_source.n) annotation (points=[-70,30; -70,40],
            style(color=3, rgbcolor={0,0,255}));
        connect(x_set.y, feed_x.s_ref) annotation (points=[59,0; 20,0; 20,50;
              12,50],            style(color=74, rgbcolor={0,0,127}));
        connect(simpleMagnet.p, u_source1.p) annotation (points=[-40,-44; -50,
              -44; -50,-40; -70,-40], style(color=3, rgbcolor={0,0,255}));
        connect(simpleMagnet.n, u_source1.n) annotation (points=[-40,-56; -50,
              -56; -50,-60; -70,-60], style(color=3, rgbcolor={0,0,255}));
        connect(simpleMagnet.flange, feed_x1.flange_b)       annotation (points=[-20,-50;
              -10,-50],         style(color=58, rgbcolor={0,127,0}));
        connect(advancedMagnet.n, u_source.n) annotation (points=[-40,44; -50,
              44; -50,40; -70,40],   style(color=3, rgbcolor={0,0,255}));
        connect(feed_x1.s_ref,x_set. y)       annotation (points=[12,-50; 20,
              -50; 20,0; 59,0], style(color=74, rgbcolor={0,0,127}));
        connect(x_set.y,comparisonWithFEA. u) annotation (points=[59,0; 50,0;
              50,70; 58,70], style(color=74, rgbcolor={0,0,127}));
        connect(feed_x.flange_b, advancedMagnet.flange) annotation (points=[-10,
              50; -20,50], style(color=58, rgbcolor={0,127,0}));
        connect(u_source.p, advancedMagnet.p) annotation (points=[-70,60; -50,
              60; -50,56; -40,56], style(color=3, rgbcolor={0,0,255}));
        connect(ground1.p, u_source1.n) annotation (points=[-70,-70; -70,-60],
            style(color=3, rgbcolor={0,0,255}));
      end ComparisonQuasiStationary;

      model ComparisonPullInStroke
      "Pull-in stroke of both solenoid models after a voltage step at time t=0"

        extends Modelica.Icons.Example;

        parameter SI.Voltage u_step = 12 "Applied voltage";

        Modelica.Electrical.Analog.Sources.StepVoltage u_source(V=u_step)
          annotation (extent=[-70,20; -50,40], rotation=270);
        annotation (Diagram,                                       experiment(StopTime=
                0.05, Tolerance=1e-007),
          Documentation(info="<html>
<p>
Have a look at <a href=\"ElectromagneticActuator\">ElectromagneticActuator</a> for general comments and at <a href=\"SimpleSolenoidModel\">SimpleSolenoidModel</a> and <a href=\"AdvancedSolenoidModel\">AdvancedSolenoidModel</a> for a detailed description of both magnetic network models.
</p>
<p>
A voltage step is applied to both solenoid models at time t=0. The armatures of both models and therewith connected loads are pulled from their rest position at maximum air gap length to their minimum position that is due to a stopper. As a reference, simulation results obtained with a dynamic model based on stationary FEA look-up tables (not part of this library) are included. Note that these reference results are valid for the default supply voltage u_step=12V DC and the default load mass m_load=0.01kg only!
</p>
<p>
Set the <b>tolerance</b> to <b>1e-7</b> and <b>simulate for 0.05 s</b>. Plot actuator current, force and position of the two magnetic network models and the FEA-based reference <b>vs. time</b> (each quantity in one common plot window):</p>
<pre>
Plot window for current:
    simpleMagnet.p.i            // rapid current rise indicates low inductance of simple network model
    advancedMagnet.p.i          // current rise slower, better match with FEA reference
    comparisonWithFEA.y[1]      // current obtained from dynamic model based on stationary FEA look-up tables

Plot window for force:
    simpleMagnet.armature.flange_a.f       // reluctance force of simple actuator model
    advancedMagnet.armature.flange_a.f     // reluctance force of advanced actuator model
    comparisonWithFEA.y[2]      // force obtained from dynamic model based on stationary FEA look-up tables

Plot window for position:
    simpleMagnet.x              // armature position of simple actuator model
    advancedMagnet.x            // armature position of advanced actuator model
    comparisonWithFEA.y[3]      // position obtained from dynamic model based on stationary FEA look-up tables
</pre>
<p>
The characteristic current drop during pull-in is due to both armature motion and increasing inductance with decreasing air gap length. Bouncing occurs when  armature and load of each model arrive at the stopper at minimum position. Although the pull-in times of the two magnetic network models are relatively close to the time obtained with the reference model, the accuracy of the advanced solenoid model is better, as one can tell from a comparison of the current rise at the beginning of the stroke.
</p>

</html>"));
        Modelica.Electrical.Analog.Basic.Ground ground
          annotation (extent=[-70,-10; -50,10]);
        Modelica.Mechanics.Translational.SlidingMass m_load(m=0.01)
        "translatory load to be pulled horizontally"
          annotation (extent=[20,20; 40,40]);
        Modelica.Blocks.Sources.CombiTimeTable comparisonWithFEA(
          table=[0,0,0,0.005; 2.61165e-007,7.93537e-005,-1.97914e-005,0.005;
              2.61165e-007,7.93537e-005,-1.97914e-005,0.005; 0.0001,0.0300045,-0.00748335,
              0.005; 0.0002,0.05926,-0.0147799,0.005; 0.0003,0.0877841,-0.021894,
              0.00499999; 0.0004,0.115593,-0.036608,0.00499997; 0.0005,0.142707,
              -0.0568957,0.00499994; 0.0006,0.169143,-0.076676,0.00499988;
              0.0007,0.194915,-0.0959614,0.0049998; 0.0008,0.220042,-0.124763,
              0.00499968; 0.0009,0.244539,-0.155317,0.00499951; 0.001,0.26842,-0.185107,
              0.00499928; 0.0011,0.291701,-0.214153,0.00499898; 0.0012,0.314394,
              -0.249655,0.0049986; 0.0013,0.336514,-0.288306,0.00499812; 0.0014,
              0.358074,-0.325991,0.00499754; 0.0015,0.379086,-0.362735,
              0.00499682; 0.0016,0.399562,-0.398563,0.00499597; 0.0017,0.419514,
              -0.44324,0.00499496; 0.0018,0.438955,-0.487015,0.00499378; 0.0019,
              0.457893,-0.529698,0.00499242; 0.002,0.47634,-0.571317,0.00499085;
              0.0021,0.494305,-0.611901,0.00498906; 0.0022,0.511799,-0.657374,
              0.00498704; 0.0023,0.528832,-0.704491,0.00498476; 0.0024,0.545412,
              -0.750434,0.00498221; 0.0025,0.561548,-0.795237,0.00497937;
              0.0026,0.577248,-0.83893,0.00497623; 0.0027,0.592521,-0.881543,
              0.00497277; 0.0028,0.607375,-0.926803,0.00496896; 0.0029,0.62182,
              -0.974598,0.0049648; 0.003,0.63586,-1.02121,0.00496027; 0.0031,
              0.649503,-1.06667,0.00495534; 0.0032,0.662756,-1.11102,0.00495;
              0.0033,0.675625,-1.15428,0.00494424; 0.0034,0.688119,-1.19648,
              0.00493803; 0.0035,0.700242,-1.23778,0.00493136; 0.0036,0.712005,
              -1.28391,0.00492421; 0.0037,0.72341,-1.32891,0.00491657; 0.0038,
              0.734463,-1.3728,0.00490842; 0.0039,0.74517,-1.41563,0.00489974;
              0.004,0.755536,-1.45743,0.00489052; 0.0041,0.765568,-1.49822,
              0.00488074; 0.0042,0.775269,-1.53803,0.00487038; 0.0043,0.784646,
              -1.57689,0.00485943; 0.0044,0.793704,-1.61483,0.00484787; 0.0045,
              0.80245,-1.65314,0.00483569; 0.0046,0.810888,-1.69366,0.00482288;
              0.0047,0.81902,-1.7332,0.00480941; 0.0048,0.826851,-1.77179,
              0.00479528; 0.0049,0.834387,-1.80945,0.00478046; 0.005,0.841631,-1.84622,
              0.00476495; 0.0051,0.84859,-1.88259,0.00474873; 0.0052,0.855304,-1.92429,
              0.00473179; 0.0053,0.861739,-1.96564,0.0047141; 0.0054,0.8679,-2.00668,
              0.00469566; 0.0055,0.873791,-2.04743,0.00467645; 0.0056,0.879419,
              -2.08794,0.00465645; 0.0057,0.884782,-2.1282,0.00463565; 0.0058,
              0.889885,-2.16824,0.00461403; 0.0059,0.894731,-2.20808,0.00459157;
              0.006,0.899322,-2.24774,0.00456827; 0.0061,0.903661,-2.28927,
              0.0045441; 0.0062,0.907752,-2.33091,0.00451905; 0.0063,0.911603,-2.37014,
              0.0044931; 0.0064,0.915232,-2.40274,0.00446624; 0.0065,0.91862,-2.43469,
              0.00443846; 0.0066,0.92177,-2.466,0.00440974; 0.0067,0.924686,-2.49668,
              0.00438007; 0.0068,0.927368,-2.52672,0.00434945; 0.0069,0.929822,
              -2.55615,0.00431785; 0.007,0.93205,-2.58498,0.00428527; 0.0071,
              0.934052,-2.61318,0.00425169; 0.0072,0.935241,-2.64973,0.00421711;
              0.0073,0.936164,-2.68643,0.00418151; 0.0074,0.936854,-2.7228,
              0.00414488; 0.0075,0.937309,-2.7588,0.0041072; 0.0076,0.937532,-2.7944,
              0.00406845; 0.0077,0.937522,-2.82958,0.00402864; 0.0078,0.937411,
              -2.866,0.00398773; 0.0079,0.937385,-2.90613,0.00394572; 0.008,
              0.937133,-2.94589,0.0039026; 0.0081,0.936656,-2.98525,0.00385834;
              0.0082,0.935953,-3.02414,0.00381293; 0.0083,0.935024,-3.06251,
              0.00376636; 0.0084,0.934308,-3.10824,0.00371862; 0.0085,0.933608,
              -3.15783,0.00366967; 0.0086,0.93269,-3.20708,0.00361952; 0.0087,
              0.931553,-3.25592,0.00356812; 0.0088,0.930194,-3.30427,0.00351548;
              0.0089,0.928473,-3.35247,0.00346157; 0.009,0.926467,-3.40014,
              0.00340636; 0.0091,0.924232,-3.44698,0.00334985; 0.0092,0.921766,
              -3.49289,0.00329202; 0.0093,0.918579,-3.53879,0.00323283; 0.0094,
              0.913925,-3.5856,0.00317229; 0.0095,0.909004,-3.63034,0.00311037;
              0.0096,0.903809,-3.67275,0.00304706; 0.0097,0.89859,-3.72881,
              0.00298233; 0.0098,0.893783,-3.82589,0.00291616; 0.0099,0.888707,
              -3.92096,0.00284852; 0.01,0.883343,-4.01357,0.00277938; 0.0101,
              0.876979,-4.10734,0.00270869; 0.0102,0.869783,-4.19987,0.00263642;
              0.0103,0.862246,-4.28752,0.00256254; 0.0104,0.854574,-4.37627,
              0.00248701; 0.0105,0.847614,-4.49154,0.00240979; 0.0106,0.840302,
              -4.60102,0.00233085; 0.0107,0.832625,-4.70399,0.00225014; 0.0108,
              0.822938,-4.82647,0.00216761; 0.0109,0.812813,-4.93752,0.00208323;
              0.011,0.802204,-5.04175,0.00199695; 0.0111,0.78997,-5.30274,
              0.00190873; 0.0112,0.777197,-5.54515,0.00181846; 0.0113,0.763521,
              -5.78149,0.00172606; 0.0114,0.748272,-6.039,0.00163144; 0.0115,
              0.73235,-6.25778,0.0015345; 0.0116,0.715211,-6.57852,0.00143514;
              0.0117,0.696998,-6.91971,0.00133326; 0.0118,0.677065,-7.30735,
              0.00122872; 0.0119,0.652791,-7.88085,0.00112136; 0.012,0.62734,-8.29718,
              0.00101097; 0.0121,0.597125,-9.13179,0.000897364; 0.0122,0.564919,
              -9.82427,0.000780251; 0.0123,0.527838,-11.1684,0.000659331;
              0.0124,0.487477,-12.1609,0.000534142; 0.0125,0.436631,-14.9103,
              0.000404205; 0.0126,0.379243,-16.2449,0.000268616; 0.0126134,
              0.371242,-16.2777,0.00025; 0.0126134,0.371242,-16.2777,0.00025;
              0.0126868,0.350822,-16.2554,0.000198624; 0.0126868,0.350822,-16.2554,
              0.000198624; 0.0127,0.351869,-16.3218,0.000199455; 0.0128,0.37695,
              -17.0338,0.000241587; 0.0128157,0.381787,-17.1198,0.00025;
              0.0128157,0.381787,-17.1198,0.00025; 0.0129,0.406591,-17.48,
              0.000292352; 0.013,0.433421,-17.8191,0.000336402; 0.0131,0.457261,
              -17.8337,0.000373609; 0.0132,0.477911,-17.6706,0.000403962;
              0.0133,0.495294,-17.4605,0.00042752; 0.0134,0.509353,-17.3988,
              0.000444358; 0.0135,0.520015,-17.4878,0.0004545; 0.0136,0.527192,
              -17.7433,0.000457911; 0.0136003,0.527207,-17.7443,0.000457911;
              0.0136003,0.527207,-17.7443,0.000457911; 0.0137,0.530748,-18.1997,
              0.000454491; 0.0138,0.530517,-18.8646,0.000444064; 0.0139,
              0.526294,-19.7142,0.000426376; 0.014,0.517828,-20.6871,
              0.000401101; 0.0141,0.504836,-21.6765,0.000367869; 0.0142,
              0.487037,-22.6627,0.000326301; 0.0143,0.464073,-23.4017,
              0.000276025; 0.0143458,0.451744,-23.5657,0.00025; 0.0143458,
              0.451744,-23.5657,0.00025; 0.0144,0.439383,-23.6302,0.000223375;
              0.0144518,0.438001,-23.8106,0.00021654; 0.0144518,0.438001,-23.8106,
              0.00021654; 0.0145,0.442437,-24.0882,0.000220288; 0.0146,0.459291,
              -24.7355,0.000241352; 0.014643,0.466338,-24.9736,0.00025;
              0.014643,0.466338,-24.9736,0.00025; 0.0147,0.47417,-25.2545,
              0.000258795; 0.0148,0.483493,-25.7045,0.000266567; 0.0148288,
              0.485111,-25.8323,0.00026698; 0.0148288,0.485111,-25.8323,
              0.00026698; 0.0149,0.486998,-26.1506,0.000264454; 0.015,0.484444,
              -26.5924,0.000252282; 0.0150127,0.483671,-26.6456,0.00025;
              0.0150127,0.483671,-26.6456,0.00025; 0.0151,0.477935,-26.9803,
              0.000233764; 0.0151954,0.478678,-27.3825,0.000227777; 0.0151954,
              0.478678,-27.3825,0.000227777; 0.0152,0.478896,-27.404,
              0.000227786; 0.0153,0.486112,-27.9096,0.000231723; 0.0154,
              0.494618,-28.4114,0.000237745; 0.0154716,0.499054,-28.7526,
              0.000239402; 0.0154716,0.499054,-28.7526,0.000239402; 0.0155,
              0.500242,-28.8872,0.000239151; 0.0156,0.502893,-29.3755,
              0.000235871; 0.0157,0.505639,-29.8643,0.000232816; 0.0158,
              0.509736,-30.3772,0.000231912; 0.0158118,0.51029,-30.4396,
              0.000231905; 0.0158118,0.51029,-30.4396,0.000231905; 0.0159,
              0.514622,-30.9065,0.000232198; 0.016,0.519654,-31.4343,
              0.000232755; 0.016048,0.521947,-31.6846,0.000232849; 0.016048,
              0.521947,-31.6846,0.000232849; 0.0161,0.524291,-31.9527,
              0.000232753; 0.0162,0.528618,-32.4638,0.000232328; 0.0163,0.53296,
              -32.9726,0.000231976; 0.0164,0.537374,-33.4793,0.000231787;
              0.0165,0.541801,-33.9827,0.000231672; 0.0166,0.546199,-34.4828,
              0.000231561; 0.0167,0.550555,-34.9795,0.000231435; 0.0168,
              0.554875,-35.4729,0.0002313; 0.0169,0.559164,-35.9631,0.000231166;
              0.017,0.56344,-36.4518,0.000231035; 0.0171,0.567726,-36.9417,
              0.000230906; 0.0172,0.571982,-37.4284,0.000230779; 0.0173,
              0.576209,-37.9119,0.000230653; 0.0174,0.580407,-38.3923,
              0.000230528; 0.0175,0.584575,-38.8695,0.000230405; 0.0176,
              0.588716,-39.3436,0.000230284; 0.0177,0.593137,-39.8493,
              0.000230163; 0.0178,0.59757,-40.357,0.000230038; 0.0179,0.601967,
              -40.8716,0.000229911; 0.018,0.60633,-41.3953,0.000229783; 0.0181,
              0.610659,-41.9153,0.000229654; 0.0182,0.614955,-42.4317,
              0.000229526; 0.0183,0.619218,-42.9441,0.0002294; 0.0184,0.623441,
              -43.452,0.000229276; 0.0185,0.627634,-43.9562,0.000229154; 0.0186,
              0.631795,-44.4569,0.000229034; 0.0187,0.635926,-44.954,
              0.000228915; 0.0188,0.640026,-45.4476,0.000228797; 0.0189,
              0.644096,-45.9377,0.000228681; 0.019,0.648136,-46.4242,
              0.000228566; 0.0191,0.652146,-46.9074,0.000228453; 0.0192,
              0.656126,-47.387,0.000228341; 0.0193,0.660077,-47.8633,
              0.000228231; 0.0194,0.663999,-48.3362,0.000228122; 0.0195,
              0.667892,-48.8057,0.000228014; 0.0196,0.671756,-49.2718,
              0.000227908; 0.0197,0.675592,-49.7347,0.000227802; 0.0198,0.67979,
              -50.2404,0.000227697; 0.0199,0.684118,-50.7623,0.000227586; 0.02,
              0.688404,-51.2799,0.000227471; 0.0201,0.692654,-51.7933,
              0.000227355; 0.0202,0.696868,-52.3025,0.000227241; 0.0203,
              0.701047,-52.8002,0.00022713; 0.0204,0.705193,-53.2717,
              0.000227022; 0.0205,0.709307,-53.7394,0.000226918; 0.0206,
              0.713479,-54.2135,0.000226817; 0.0207,0.717635,-54.686,
              0.000226716; 0.0208,0.721755,-55.1544,0.000226615; 0.0209,
              0.725839,-55.619,0.000226515; 0.021,0.729888,-56.0796,0.000226416;
              0.0211,0.733903,-56.5364,0.000226319; 0.0212,0.737883,-56.9893,
              0.000226222; 0.0213,0.741829,-57.4383,0.000226127; 0.0214,
              0.745732,-57.8827,0.000226033; 0.0215,0.749587,-58.3217,
              0.000225941; 0.0216,0.75341,-58.7569,0.00022585; 0.0217,0.757199,
              -59.1885,0.00022576; 0.0218,0.760956,-59.6164,0.000225671; 0.0219,
              0.764681,-60.0407,0.000225583; 0.022,0.768373,-60.4614,
              0.000225497; 0.0221,0.772034,-60.8786,0.000225411; 0.0222,
              0.775663,-61.2922,0.000225326; 0.0223,0.779579,-61.7378,
              0.000225242; 0.0224,0.784355,-62.2802,0.000225151; 0.0225,
              0.789065,-62.8168,0.000225046; 0.0226,0.793716,-63.3474,
              0.000224938; 0.0227,0.798315,-63.8721,0.000224831; 0.0228,
              0.802863,-64.3256,0.000224728; 0.0229,0.80737,-64.7356,
              0.000224637; 0.023,0.811833,-65.1406,0.000224555; 0.0231,0.816247,
              -65.541,0.000224477; 0.0232,0.820611,-65.9369,0.000224399; 0.0233,
              0.824909,-66.3269,0.000224322; 0.0234,0.829106,-66.7079,
              0.000224246; 0.0235,0.833258,-67.0845,0.000224172; 0.0236,
              0.837362,-67.457,0.000224099; 0.0237,0.84142,-67.8252,0.000224027;
              0.0238,0.845433,-68.1893,0.000223957; 0.0239,0.8494,-68.5494,
              0.000223887; 0.024,0.853323,-68.9053,0.000223818; 0.0241,0.857201,
              -69.2573,0.00022375; 0.0242,0.861036,-69.6053,0.000223683; 0.0243,
              0.864828,-69.9494,0.000223617; 0.0244,0.868577,-70.2896,
              0.000223552; 0.0245,0.873541,-70.7381,0.000223484; 0.0246,
              0.878506,-71.1879,0.000223404; 0.0247,0.883389,-71.6312,
              0.00022332; 0.0248,0.888198,-72.0678,0.000223236; 0.0249,0.892935,
              -72.4978,0.000223154; 0.025,0.8976,-72.9212,0.000223074; 0.0251,
              0.902194,-73.2832,0.000222997; 0.0252,0.906729,-73.5797,
              0.00022293; 0.0253,0.911203,-73.8711,0.000222873; 0.0254,0.915611,
              -74.1579,0.00022282; 0.0255,0.919953,-74.4404,0.000222768; 0.0256,
              0.924227,-74.7185,0.000222716; 0.0257,0.928436,-74.9923,
              0.000222665; 0.0258,0.932872,-75.2805,0.000222615; 0.0259,
              0.937419,-75.5759,0.000222563; 0.026,0.941886,-75.8664,
              0.000222509; 0.0261,0.946276,-76.1519,0.000222456; 0.0262,
              0.950592,-76.4326,0.000222404; 0.0263,0.954834,-76.7084,
              0.000222354; 0.0264,0.959005,-76.9795,0.000222304; 0.0265,
              0.963104,-77.246,0.000222255; 0.0266,0.967134,-77.5079,
              0.000222207; 0.0267,0.971094,-77.7654,0.00022216; 0.0268,0.974988,
              -78.0184,0.000222114; 0.0269,0.978815,-78.2671,0.000222068; 0.027,
              0.982577,-78.5115,0.000222024; 0.0271,0.986275,-78.7518,
              0.00022198; 0.0272,0.98991,-78.9879,0.000221937; 0.0273,0.993484,
              -79.2201,0.000221895; 0.0274,0.996996,-79.4482,0.000221854;
              0.0275,1.00082,-79.6845,0.000221813; 0.0276,1.00486,-79.8903,
              0.000221773; 0.0277,1.00883,-80.0919,0.000221735; 0.0278,1.01272,
              -80.2892,0.000221699; 0.0279,1.01653,-80.4824,0.000221665; 0.028,
              1.02026,-80.6717,0.000221631; 0.0281,1.02392,-80.8572,0.000221597;
              0.0282,1.0275,-81.0389,0.000221565; 0.0283,1.03101,-81.2168,
              0.000221533; 0.0284,1.03445,-81.3911,0.000221501; 0.0285,1.03781,
              -81.5619,0.000221471; 0.0286,1.04111,-81.7292,0.000221441; 0.0287,
              1.04434,-81.893,0.000221412; 0.0288,1.04751,-82.0535,0.000221383;
              0.0289,1.05061,-82.2107,0.000221355; 0.029,1.05365,-82.3647,
              0.000221328; 0.0291,1.05663,-82.5155,0.000221301; 0.0292,1.05954,
              -82.6633,0.000221275; 0.0293,1.0624,-82.808,0.000221249; 0.0294,
              1.0652,-82.9498,0.000221224; 0.0295,1.06794,-83.0887,0.000221199;
              0.0296,1.07063,-83.2248,0.000221175; 0.0297,1.07326,-83.3581,
              0.000221151; 0.0298,1.07584,-83.4886,0.000221128; 0.0299,1.07836,
              -83.6165,0.000221106; 0.03,1.08088,-83.7439,0.000221083; 0.0301,
              1.08376,-83.8895,0.000221061; 0.0302,1.08657,-84.0316,0.000221037;
              0.0303,1.08931,-84.1703,0.000221012; 0.0304,1.09198,-84.3057,
              0.000220988; 0.0305,1.09459,-84.4378,0.000220965; 0.0306,1.09714,
              -84.5667,0.000220942; 0.0307,1.09962,-84.6924,0.00022092; 0.0308,
              1.10205,-84.7987,0.000220899; 0.0309,1.10442,-84.8994,0.00022088;
              0.031,1.10673,-84.9975,0.000220862; 0.0311,1.10898,-85.0932,
              0.000220846; 0.0312,1.11119,-85.1866,0.000220829; 0.0313,1.11333,
              -85.2778,0.000220813; 0.0314,1.11543,-85.3668,0.000220798; 0.0315,
              1.11748,-85.4536,0.000220782; 0.0316,1.11947,-85.5382,0.000220768;
              0.0317,1.12142,-85.6209,0.000220753; 0.0318,1.12332,-85.7015,
              0.000220739; 0.0319,1.12518,-85.7802,0.000220725; 0.032,1.12699,-85.857,
              0.000220712; 0.0321,1.12875,-85.9319,0.000220699; 0.0322,1.13048,
              -86.005,0.000220686; 0.0323,1.13216,-86.0763,0.000220673; 0.0324,
              1.1338,-86.1459,0.000220661; 0.0325,1.1354,-86.2138,0.000220649;
              0.0326,1.13696,-86.28,0.000220638; 0.0327,1.13849,-86.3447,
              0.000220627; 0.0328,1.13997,-86.4078,0.000220616; 0.0329,1.14143,
              -86.4693,0.000220605; 0.033,1.14284,-86.5294,0.000220594; 0.0331,
              1.14423,-86.588,0.000220584; 0.0332,1.14558,-86.6452,0.000220574;
              0.0333,1.14689,-86.701,0.000220564; 0.0334,1.14818,-86.7555,
              0.000220555; 0.0335,1.14943,-86.8086,0.000220546; 0.0336,1.15065,
              -86.8605,0.000220537; 0.0337,1.15185,-86.9111,0.000220528; 0.0338,
              1.15301,-86.9605,0.000220519; 0.0339,1.15415,-87.0086,0.000220511;
              0.034,1.15526,-87.0556,0.000220503; 0.0341,1.15634,-87.1015,
              0.000220495; 0.0342,1.1574,-87.1463,0.000220487; 0.0343,1.15843,-87.19,
              0.000220479; 0.0344,1.15943,-87.2326,0.000220472; 0.0345,1.16041,
              -87.2742,0.000220465; 0.0346,1.16137,-87.3148,0.000220458; 0.0347,
              1.16231,-87.3544,0.000220451; 0.0348,1.16322,-87.3931,0.000220444;
              0.0349,1.16411,-87.4308,0.000220438; 0.035,1.16498,-87.4676,
              0.000220431; 0.0351,1.16582,-87.5035,0.000220425; 0.0352,1.16665,
              -87.5385,0.000220419; 0.0353,1.16746,-87.5727,0.000220413; 0.0354,
              1.16824,-87.6061,0.000220407; 0.0355,1.16901,-87.6386,0.000220402;
              0.0356,1.16976,-87.6704,0.000220396; 0.0357,1.17049,-87.7014,
              0.000220391; 0.0358,1.17121,-87.7316,0.000220386; 0.0359,1.1719,-87.7612,
              0.00022038; 0.036,1.17258,-87.79,0.000220375; 0.0361,1.17325,-87.8181,
              0.000220371; 0.0362,1.1739,-87.8455,0.000220366; 0.0363,1.17453,-87.8722,
              0.000220361; 0.0364,1.17514,-87.8984,0.000220357; 0.0365,1.17574,
              -87.9238,0.000220352; 0.0366,1.17633,-87.9487,0.000220348; 0.0367,
              1.1769,-87.9729,0.000220344; 0.0368,1.17746,-87.9966,0.00022034;
              0.0369,1.17801,-88.0197,0.000220336; 0.037,1.17858,-88.0441,
              0.000220332; 0.0371,1.17922,-88.0712,0.000220328; 0.0372,1.17985,
              -88.0975,0.000220323; 0.0373,1.18045,-88.123,0.000220319; 0.0374,
              1.18103,-88.1477,0.000220314; 0.0375,1.1816,-88.1717,0.00022031;
              0.0376,1.18215,-88.195,0.000220306; 0.0377,1.18268,-88.2176,
              0.000220302; 0.0378,1.1832,-88.2395,0.000220299; 0.0379,1.1837,-88.2607,
              0.000220295; 0.038,1.18419,-88.2814,0.000220291; 0.0381,1.18466,-88.3014,
              0.000220288; 0.0382,1.18512,-88.3208,0.000220284; 0.0383,1.18556,
              -88.3396,0.000220281; 0.0384,1.18599,-88.3578,0.000220278; 0.0385,
              1.18641,-88.3756,0.000220275; 0.0386,1.18682,-88.3928,0.000220272;
              0.0387,1.18721,-88.4094,0.000220269; 0.0388,1.18759,-88.4256,
              0.000220266; 0.0389,1.18796,-88.4413,0.000220264; 0.039,1.18832,-88.4565,
              0.000220261; 0.0391,1.18867,-88.4713,0.000220258; 0.0392,1.18901,
              -88.4856,0.000220256; 0.0393,1.18934,-88.4995,0.000220253; 0.0394,
              1.18965,-88.513,0.000220251; 0.0395,1.18996,-88.5261,0.000220249;
              0.0396,1.19026,-88.5388,0.000220247; 0.0397,1.19055,-88.5511,
              0.000220245; 0.0398,1.19084,-88.563,0.000220242; 0.0399,1.19111,-88.5746,
              0.00022024; 0.04,1.19137,-88.5859,0.000220239; 0.0401,1.19163,-88.5968,
              0.000220237; 0.0402,1.19188,-88.6074,0.000220235; 0.0403,1.19212,
              -88.6176,0.000220233; 0.0404,1.19236,-88.6276,0.000220231; 0.0405,
              1.19259,-88.6373,0.00022023; 0.0406,1.19281,-88.6466,0.000220228;
              0.0407,1.19302,-88.6557,0.000220226; 0.0408,1.19323,-88.6646,
              0.000220225; 0.0409,1.19343,-88.6731,0.000220223; 0.041,1.19363,-88.6814,
              0.000220222; 0.0411,1.19382,-88.6895,0.000220221; 0.0412,1.19401,
              -88.6973,0.000220219; 0.0413,1.19418,-88.7049,0.000220218; 0.0414,
              1.19436,-88.7122,0.000220217; 0.0415,1.19453,-88.7194,0.000220215;
              0.0416,1.19469,-88.7263,0.000220214; 0.0417,1.19485,-88.733,
              0.000220213; 0.0418,1.195,-88.7395,0.000220212; 0.0419,1.19515,-88.7459,
              0.000220211; 0.042,1.1953,-88.752,0.00022021; 0.0421,1.19544,-88.7579,
              0.000220209; 0.0422,1.19557,-88.7637,0.000220208; 0.0423,1.19571,
              -88.7693,0.000220207; 0.0424,1.19583,-88.7747,0.000220206; 0.0425,
              1.19596,-88.78,0.000220205; 0.0426,1.19608,-88.7851,0.000220204;
              0.0427,1.1962,-88.7901,0.000220203; 0.0428,1.19631,-88.7949,
              0.000220202; 0.0429,1.19642,-88.7996,0.000220202; 0.043,1.19653,-88.8041,
              0.000220201; 0.0431,1.19663,-88.8085,0.0002202; 0.0432,1.19673,-88.8127,
              0.000220199; 0.0433,1.19683,-88.8169,0.000220199; 0.0434,1.19692,
              -88.8209,0.000220198; 0.0435,1.19702,-88.8248,0.000220197; 0.0436,
              1.1971,-88.8286,0.000220197; 0.0437,1.19719,-88.8322,0.000220196;
              0.0438,1.19728,-88.8358,0.000220195; 0.0439,1.19736,-88.8392,
              0.000220195; 0.044,1.19744,-88.8426,0.000220194; 0.0441,1.19751,-88.8458,
              0.000220194; 0.0442,1.19759,-88.8489,0.000220193; 0.0443,1.19766,
              -88.852,0.000220192; 0.0444,1.19773,-88.855,0.000220192; 0.0445,
              1.1978,-88.8578,0.000220191; 0.0446,1.19786,-88.8606,0.000220191;
              0.0447,1.19793,-88.8633,0.000220191; 0.0448,1.19799,-88.8659,
              0.00022019; 0.0449,1.19805,-88.8685,0.00022019; 0.045,1.19811,-88.871,
              0.000220189; 0.0451,1.19816,-88.8734,0.000220189; 0.0452,1.19822,
              -88.8757,0.000220188; 0.0453,1.19827,-88.8779,0.000220188; 0.0454,
              1.19832,-88.8801,0.000220188; 0.0455,1.19837,-88.8822,0.000220187;
              0.0456,1.19842,-88.8843,0.000220187; 0.0457,1.19847,-88.8863,
              0.000220187; 0.0458,1.19851,-88.8882,0.000220186; 0.0459,1.19856,
              -88.8901,0.000220186; 0.046,1.1986,-88.8919,0.000220186; 0.0461,
              1.19864,-88.8937,0.000220185; 0.0462,1.19868,-88.8954,0.000220185;
              0.0463,1.19872,-88.897,0.000220185; 0.0464,1.19876,-88.8987,
              0.000220184; 0.0465,1.1988,-88.9002,0.000220184; 0.0466,1.19883,-88.9017,
              0.000220184; 0.0467,1.19887,-88.9032,0.000220184; 0.0468,1.1989,-88.9046,
              0.000220183; 0.0469,1.19893,-88.906,0.000220183; 0.047,1.19897,-88.9074,
              0.000220183; 0.0471,1.199,-88.9087,0.000220183; 0.0472,1.19903,-88.91,
              0.000220182; 0.0473,1.19906,-88.9112,0.000220182; 0.0474,1.19908,
              -88.9124,0.000220182; 0.0475,1.19911,-88.9135,0.000220182; 0.0476,
              1.19914,-88.9146,0.000220182; 0.0477,1.19916,-88.9157,0.000220181;
              0.0478,1.19919,-88.9168,0.000220181; 0.0479,1.19921,-88.9178,
              0.000220181; 0.048,1.19924,-88.9188,0.000220181; 0.0481,1.19926,-88.9198,
              0.000220181; 0.0482,1.19928,-88.9207,0.000220181; 0.0483,1.1993,-88.9216,
              0.00022018; 0.0484,1.19932,-88.9225,0.00022018; 0.0485,1.19934,-88.9233,
              0.00022018; 0.0486,1.19936,-88.9242,0.00022018; 0.0487,1.19938,-88.925,
              0.00022018; 0.0488,1.1994,-88.9258,0.00022018; 0.0489,1.19942,-88.9265,
              0.00022018; 0.049,1.19944,-88.9273,0.000220179; 0.0491,1.19945,-88.928,
              0.000220179; 0.0492,1.19947,-88.9287,0.000220179; 0.0493,1.19948,
              -88.9293,0.000220179; 0.0494,1.1995,-88.93,0.000220179; 0.0495,
              1.19951,-88.9306,0.000220179; 0.0496,1.19953,-88.9312,0.000220179;
              0.0497,1.19954,-88.9318,0.000220179; 0.0498,1.19956,-88.9324,
              0.000220179; 0.0499,1.19957,-88.933,0.000220178; 0.05,1.19958,-88.9335,
              0.000220178; 0.05,1.19958,-88.9335,0.000220178],
          tableOnFile=false,
          columns=2:4,
          extrapolation=Modelica.Blocks.Types.Extrapolation.HoldLastPoint)
        "Valid for u_source=12VDC and m_load=0.01kg only; column 2: current, col.3: force, col.4: position"
          annotation (extent=[-40,60; -20,80]);
        Modelica_Magnetic.Examples.ElectromagneticActuator.AdvancedSolenoidModel
        advancedMagnet               annotation (extent=[-20,20; 0,40]);
        Modelica.Electrical.Analog.Sources.StepVoltage u_source1(V=u_step)
          annotation (extent=[-70,-60; -50,-40],
                                               rotation=270);
        Modelica.Electrical.Analog.Basic.Ground ground1
          annotation (extent=[-70,-90; -50,-70]);
        Modelica.Mechanics.Translational.SlidingMass m_load1(
                                                            m=0.01)
        "translatory load to be pulled horizontally"
          annotation (extent=[20,-60; 40,-40]);
        Modelica_Magnetic.Examples.ElectromagneticActuator.SimpleSolenoidModel
        simpleMagnet                 annotation (extent=[-20,-60; 0,-40]);
      equation
        connect(ground.p, u_source.n) annotation (points=[-60,10; -60,20],
            style(color=3, rgbcolor={0,0,255}));
        connect(u_source.p,advancedMagnet. p)
                                      annotation (points=[-60,40; -40,40; -40,
              36; -20,36],
                       style(color=3, rgbcolor={0,0,255}));
        connect(advancedMagnet.n, u_source.n)
                                      annotation (points=[-20,24; -40,24; -40,
              20; -60,20],   style(color=3, rgbcolor={0,0,255}));
        connect(advancedMagnet.flange, m_load.flange_a)
                                                annotation (points=[0,30; 20,30],
                    style(color=58, rgbcolor={0,127,0}));
        connect(ground1.p, u_source1.n)
                                      annotation (points=[-60,-70; -60,-60],
            style(color=3, rgbcolor={0,0,255}));
        connect(u_source1.p,simpleMagnet. p)
                                      annotation (points=[-60,-40; -40,-40; -40,
              -44; -20,-44],
                       style(color=3, rgbcolor={0,0,255}));
        connect(simpleMagnet.n, u_source1.n)
                                      annotation (points=[-20,-56; -40,-56; -40,
              -60; -60,-60], style(color=3, rgbcolor={0,0,255}));
        connect(simpleMagnet.flange, m_load1.flange_a)
                                                annotation (points=[0,-50; 20,
              -50], style(color=58, rgbcolor={0,127,0}));
      end ComparisonPullInStroke;

      annotation (Documentation(info="<html>
<p>
In electromagnetic or reluctance actuators, a thrust or reluctance force is generated due to a non-zero gradient of the relative magnetic permeability my_r at surfaces between regions of different permeability (non-saturated ferromagnetic material: my_r>>1, adjacent air: my_r=1). In lumped magnetic network models, this force can be calculated as shortly outlined in <a href=\"UsersGuide.ReluctanceForceCalculation\">Reluctance Force</a> of the Users Guide.
</p>
<p>
As an example of a reluctance actuator, a simple axisymmetric lifting magnet with planar end planes of armature and pole is shown. Often, a <a href=\"SimpleSolenoidModel\">SimpleSolenoidModel</a> is sufficient for initial rough design of such an actuator's magnetic subsystem. Higher accuracy can be gained from an <a href=\"AdvancedSolenoidModel\">AdvancedSolenoidModel</a> where the coil-imposed magnetomotive force is split and the leakage flux between armature and yoke is accounted for more precisely.
</p>
<p>
The differences between these two models in static behaviour can be analysed and compared to results obtained with a more accurate finite element analysis (FEA) in <a href=\"ComparisonQuasiStationary\">ComparisonQuasiStationary</a>. The resulting differences in dynamic behaviour can be analysed and compared to FEA results with simulation of a pull-in stroke in <a href=\"ComparisonPullInStroke\">ComparisonPullInStroke</a>.
</p>
</html>"), Icon(
          Rectangle(extent=[-100,-100; 80,50],   style(fillColor=7)),
          Polygon(points=[-100,50; -80,70; 100,70; 80,50; -100,50],      style(
                fillColor=7)),
          Polygon(points=[100,70; 100,-80; 80,-100; 80,50; 100,70],      style(
                fillColor=7)),
          Text(
            extent=[-100,70; 100,-130],
            style(color=3, rgbcolor={0,0,255}),
            string="3"),
          Text(
            extent=[-120,130; 120,70],
            string="%name",
            style(color=3))));

    end ElectromagneticActuator;
  end Examples;


  package Interfaces "Interfaces of magnetic network components"
    extends Modelica.Icons.Library;

    annotation (Coordsys(
        extent=[-100, -100; 100, 100],
        grid=[2, 2],
        component=[20, 20]), Window(
        x=0,
        y=0.48,
        width=0.2,
        height=0.21,
        library=1,
        autolayout=1),
    Documentation(info="<HTML>
<p>
This package contains connectors for the magnetic domain, partial models for lumped magnetic network components and a partial model with connectors for a one-phase electro-mechanical translatory actuator.
</p>

</HTML>"));

    connector MagneticPort "Generic magnetic port"
      SI.MagneticPotential V_mag "Magnetic potential at the port";
      flow SI.MagneticFlux Phi "Magnetic flux flowing into the port";

      annotation (defaultComponentName = "mag",
    Icon(Rectangle(extent=[-100,100; 100,-100], style(
              color=45,
              fillColor=45,
              fillPattern=1))),
        Diagram(                Text(
            extent=[-160,110; 40,50],
            style(color=45),
            string="%name"),  Rectangle(extent=[-40,40; 40,-40], style(
              color=45,
              fillColor=45))));
    end MagneticPort;

    connector PositiveMagneticPort "Positive magnetic port"
      SI.MagneticPotential V_mag "Magnetic potential at the port";
      flow SI.MagneticFlux Phi "Magnetic flux flowing into the port";

      annotation (defaultComponentName = "mag_p",
    Icon(Rectangle(extent=[-100,100; 100,-100], style(
              color=45,
              fillColor=45,
              fillPattern=1))),
        Diagram(                Text(
            extent=[-160,110; 40,50],
            style(color=45),
            string="%name"),  Rectangle(extent=[-40,40; 40,-40], style(
              color=45,
              fillColor=45))));

    end PositiveMagneticPort;

    connector NegativeMagneticPort "Negative magnetic port"
      SI.MagneticPotential V_mag "Magnetic potential at the port";
      flow SI.MagneticFlux Phi "Magnetic flux flowing into the port";

    annotation (defaultComponentName = "mag_n",
    Icon(Rectangle(extent=[-100,100; 100,-100], style(
              color=45,
              fillColor=7,
              fillPattern=1))),
        Diagram(                Text(
            extent=[-40,110; 160,50],
            string="%name",
            style(color=45)), Rectangle(extent=[-40,40; 40,-40], style(
              color=45,
              fillColor=7))));

    end NegativeMagneticPort;

    partial model TwoMagneticPorts
    "Partial component with two magnetic ports p and n"

      PositiveMagneticPort p "Positive magnetic port"
        annotation (extent=[-110,-10; -90,10]);
      NegativeMagneticPort n "Negative magnetic port"
        annotation (extent=[90,-10; 110,10]);

      annotation (
        Coordsys(
          extent=[-100, -100; 100, 100],
          grid=[2, 2],
          component=[20, 20]),
        Documentation(info="<HTML>
<P>
Partial model of a flux tube component with <b>two</b> magnetic ports:
the positive port connector <i>p</i>, and the negative port
connector <i>n</i>.
</P>
</HTML>
"),     Diagram,
        Window(
          x=0.33,
          y=0.04,
          width=0.63,
          height=0.67),
        Icon(Text(
            extent=[-140,90; 140,40],
            style(color=3, rgbcolor={0,0,255}),
            string="%name")));
    equation

    end TwoMagneticPorts;

    partial model TwoPortComponent
    "Partial component with magnetic potential difference between two magnetic ports p and n and magnetic flux Phi from p to n"

      extends TwoMagneticPorts;

      SI.MagneticPotentialDifference V_mag
      "Magnetic potential difference between both ports";
      SI.MagneticFlux Phi "Magnetic flux from port p to port n";

    equation
      V_mag = p.V_mag - n.V_mag;
      Phi = p.Phi;
      0 = p.Phi + n.Phi;

      annotation (Documentation(info="<html>
<p>
It is assumed that the magnetic flux flowing into pin p is identical to the flux flowing out of pin n.
This magnetic flux is provided explicitly as flux Phi.
</html>"));
    end TwoPortComponent;

    partial model ElectromechanicalActuator
    "Connectors of a generic translatory actuator with one electrical phase"

      Modelica.Electrical.Analog.Interfaces.NegativePin n
      "Electrical connector"
        annotation (extent=[-110,-70; -90,-50]);
      Modelica.Electrical.Analog.Interfaces.PositivePin p
      "Electrical connector"
        annotation (extent=[-110,50; -90,70]);
      Modelica.Mechanics.Translational.Interfaces.Flange_a flange
      "Flange for connection to translatory load"
        annotation (extent=[90,-10; 110,10]);
    annotation (Diagram, Icon(Rectangle(extent=[-90,100; 90,-100],   style(
            color=45,
            rgbcolor={255,128,0},
            fillColor=7,
            rgbfillColor={255,255,255})), Text(
          extent=[-140,150; 140,100],
          string="%name")),
        Documentation(info="<html>
<p>
This partial model defines connectors for any kind of electro-mechanical actuator with translatory motion
and one electrical phase (e.g. electrodynamic, electromagnetic, piezoelectric and thermal actuators).
</p>
</html>"));
    end ElectromechanicalActuator;
  end Interfaces;


annotation (
  Coordsys(
    extent=[-100, -100; 100, 100],
    grid=[2, 2],
    component=[20, 20]),
  Window(
    x=0.45,
    y=0.01,
    width=0.44,
    height=0.65,
    library=1,
    autolayout=1),
  Documentation(info="<html>
<p>
This library contains components for modelling of electromagnetic devices with lumped magnetic networks. Those models are suited for both rough design of the magnetic subsystem of a device as well as for efficient dynamic simulation at system level together with neighbouring subsystems. At present, components and examples for modelling of <i>translatory</i> electromagnetic and electrodynamic actuators are provided. If needed, these components can be adapted to network modellling of <i>rotational</i> electrical machines.
</p>
<p>
<a href=\"Modelica_Magnetic.UsersGuide\">Users Guide</a> gives a short introduction to the underlying concept of <b>magnetic flux tubes</b>, summarizes the calculation of magnetic <b>reluctance forces</b> from lumped magnetic network models and lists <b>reference literature</b>.
</p>
<p>
<a href=\"Modelica_Magnetic.Examples\">Examples</a> illustrates the usage of magnetic network models with simple models from different fields of application.
<br>
<br>
</p>
<p>
<dl>
<dt>
<b>Copyright:</b></dt>
<br>
<dd>
Copyright (C) 2005-2007, Modelica Association and Thomas B&ouml;drich.<br>
<br>
<i>The Modelica package is <b>free</b> software; it can be redistributed and/or modified
under the terms of the <b>Modelica license</b>, see the license conditions
and the accompanying <b>disclaimer</b> in the documentation of package Modelica in file \"Modelica/package.mo\".</i></dd>
</dl>
<br>
</HTML>
"),
  Icon(
    Rectangle(extent=[36,2; 60,-84],   style(color=0)),
    Line(points=[20, 2; -60, 2; -60, -84; 20, -84; 20, -62; -36, -62; -36, -20;
           20, -20; 20, 2], style(color=0)),
    Line(points=[-64,-28; -32,-36],   style(color=0)),
    Line(points=[-64,-52; -32,-60],   style(color=0)),
    Line(points=[-64,-16; -32,-24],   style(color=0)),
    Line(points=[-64,-16; -78,-16],   style(color=0)),
    Line(points=[-64,-52; -78,-52],   style(color=0)),
    Ellipse(extent=[-78,-20; -86,-12],   style(color=0)),
    Ellipse(extent=[-78,-56; -86,-48],   style(color=0)),
    Line(points=[-64,-40; -32,-48],   style(color=0))),
    version="1.0", versionDate="2007-10-11",
  uses(Modelica(version="2.2.1")));
  import SI = Modelica.SIunits;
  import Modelica.Constants.pi;
  import my_0 = Modelica.Constants.mue_0;


  extends Modelica.Icons.Library2;


  package FluxTube
  "Lumped elements representing flux tubes of different shapes"

    package FixedShape
    "Flux tubes with fixed shape during simulation and linear or nonlinear material characteristics"
      partial model PartialFixedShape
      "Base class for flux tubes with fixed shape during simulation; linear or nonlinear material characteristics"

        extends Modelica_Magnetic.Interfaces.TwoPortComponent;

        parameter Boolean nonLinearPermeability = true
        "= true, if non-linear rel. permeability is used, otherwise constant rel. permeability"
                                                                                                  annotation(Dialog(group="Material"));
        parameter SI.RelativePermeability my_rConst = 1
        "Constant relative permeability; used if nonLinearPermeability = false"
           annotation(Dialog(group="Material", enable = not nonLinearPermeability));

        replaceable record Material =
          Modelica_Magnetic.Material.SoftMagnetic.ApproximationData
        "Ferromagnetic material characteristics; used if nonLinearPermeability = true"
                                                  annotation(choicesAllMatching=true, Dialog(group="Material", enable = nonLinearPermeability));

        Material mat;

        SI.Reluctance R_m "Magnetic reluctance";
        SI.Permeance G_m "Magnetic permeance";
        SI.MagneticFluxDensity B "Magnetic flux density";
        SI.CrossSection A "Cross-sectional area penetrated by magnetic flux";

        annotation (
          Coordsys(
            extent=[-100, -100; 100, 100],
            grid=[2, 2],
            component=[20, 20]),
          Window(
            x=0.16,
            y=0.15,
            width=0.6,
            height=0.6),
          Icon(
            Rectangle(extent=[-70,30; 70,-30], style(
                color=45,
                rgbcolor={255,128,0},
                fillColor=7,
                rgbfillColor={255,255,255},
                fillPattern=1)),
            Line(points=[-70,0; -90,0], style(
                color=45,
                rgbcolor={255,128,0},
                fillColor=7,
                rgbfillColor={255,255,255},
                fillPattern=1)),
            Line(points=[70,0; 90,0], style(
                color=45,
                rgbcolor={255,128,0},
                fillColor=7,
                rgbfillColor={255,255,255},
                fillPattern=1))),
          Diagram,
        Documentation(info="<HTML>
<p>
Please refer to the description of  the enclosing sub-package <a href=\"FixedShape\">FixedShape</a> for a description of all elements of this package.
</p>
</HTML>"));
          SI.RelativePermeability my_r "Relative magnetic permeability";

      equation
        my_r = if nonLinearPermeability then
          Modelica_Magnetic.Material.SoftMagnetic.my_rApprox(
                B,
                mat.my_i,
                mat.B_myMax,
                mat.c_a,
                mat.c_b,
                mat.n) else my_rConst;
        R_m = 1/G_m;
        V_mag = Phi * R_m;
        B = Phi/A;

      end PartialFixedShape;

      model HollowCylinderAxialFlux
      "(Hollow) cylinder with axial flux; fixed shape; linear or nonlinear material characteristics"

        extends Modelica_Magnetic.FluxTube.FixedShape.PartialFixedShape;

        parameter SI.Length l = 0.01 "Axial length (in direction of flux)"
                                                  annotation(Dialog(group="Fixed geometry"));
        parameter SI.Radius r_i = 0
        "Inner radius of hollow cylinder (zero for cylinder)"                             annotation(Dialog(group="Fixed geometry"));
        parameter SI.Radius r_o = 0.01 "Outer radius of (hollow) cylinder" annotation(Dialog(group="Fixed geometry"));

        annotation (Images(Parameters(group="Fixed geometry", source="Images/Magnetic/FluxTube/HollowCylinderAxialFlux.png")),
          Coordsys(
            extent=[-100, -100; 100, 100],
            grid=[2, 2],
            component=[20, 20]),
          Window(
            x=0.16,
            y=0.15,
            width=0.6,
            height=0.6),
          Icon,
          Diagram,
        Documentation(info="<HTML>
<p>
Please refer to the description of  the enclosing sub-package <a href=\"FixedShape\">FixedShape</a> for a description of all elements of this package.
</p>
<p>
Set the inner radius r_i=0 for modelling of a solid cylindric flux tube.
</p>
</HTML>"));
      equation
        A = pi*(r_o^2 - r_i^2);
        G_m = (my_0 * my_r * A)/ l;

      end HollowCylinderAxialFlux;

      model HollowCylinderRadialFlux
      "Hollow cylinder with radial flux; fixed shape; linear or nonlinear material characteristics"

        extends Modelica_Magnetic.FluxTube.FixedShape.PartialFixedShape;

        parameter SI.Breadth b = 0.01 "Breadth (orthogonal to flux direction)"
                                                  annotation(Dialog(group="Fixed geometry"));
        parameter SI.Radius r_i = 0.01 "Inner radius of hollow cylinder"
                                            annotation(Dialog(group="Fixed geometry"));
        parameter SI.Radius r_o = 0.02 "Outer radius of hollow cylinder" annotation(Dialog(group="Fixed geometry"));

        annotation (Images(Parameters(group="Fixed geometry", source="Images/Magnetic/FluxTube/HollowCylinderRadialFlux_fixedShape.png")),
          Coordsys(
            extent=[-100, -100; 100, 100],
            grid=[2, 2],
            component=[20, 20]),
          Window(
            x=0.16,
            y=0.15,
            width=0.6,
            height=0.6),
          Icon,
          Diagram,
        Documentation(info="<HTML>
<p>
Please refer to the description of  the enclosing sub-package <a href=\"FixedShape\">FixedShape</a> for a description of all elements of this package.
</p>
<p>
For hollow cylindric flux tubes with a radial magnetic flux, the flux density is a function of the radius. For that reason, the characteristic my_r(B) is evaluated for the flux density at the flux tube's mean radius. </p>
<p>
For those flux tube sections of a magnetic device that have a nonlinear material characteristic my_r(B) and a large aspect ratio of outer to inner radius r_o/r_i, the section can be split up in a series connection of several hollow cylindric flux tubes with radial flux. This allows for more realistic modelling of the dependence of flux density on the radius compared to modelling with just one flux tube element.
</p>
</HTML>"));
      equation
        A = b * pi*(r_o + r_i); //Area at arithmetic mean radius for calculation of average flux density
        G_m = 2* pi* my_0* my_r* b/ Modelica.Math.log(r_o/r_i);

      end HollowCylinderRadialFlux;

      annotation (Documentation(info="<html>
<p>
Please have a look at <a href=\"UsersGuide.ReluctanceForceCalculation\">Reluctance forces</a> in the Users Guide for an explanation of the different flux tube categories and resulting sub-packages.
</p>
<p>
Due to the restrictions on reluctance force calculation outlined there, flux tube elements with a possibly nonlinear material characteristic my_r(B) must have a fixed shape during simulation of actuator motion. Hence, the dimensions of these flux tubes are defined as parameters in the model components that extend the base class <a href=\"PartialFixedShape\">PartialFixedShape</a>.  </p>
<p>
For initial design of magnetic circuits, the relative permeability of possibly nonlinear flux tube elements can easily be set to a constant value my_rConst (nonLinearPermeability set to false). In some cases, this can simplify the rough geometric design of a device's magnetic circuit. Once an initial geometry is found, the magnetic subsystem can be simulated and fine-tuned with more realistic non-linear characteristics of ferromagnetic materials. Doing so requires setting of the parameter nonLinearPermeability to true and selection of one of the soft magnetic materials of <a href=\"Material.SoftMagnetic\">Material.SoftMagnetic</a>.
</p>
</html>"));
      model Cuboid
      "Flux tube with rectangular cross-section; fixed shape; linear or nonlinear material characteristics"

        extends Modelica_Magnetic.FluxTube.FixedShape.PartialFixedShape;

        parameter SI.Length l = 0.01 "Length in direction of flux"
                                                  annotation(Dialog(group="Fixed geometry"));
        parameter SI.Breadth a = 0.01 "Breadth of rectangular cross-section"                           annotation(Dialog(group="Fixed geometry"));
        parameter SI.Height b = 0.01 "Height of rectangular cross-section" annotation(Dialog(group="Fixed geometry"));

        annotation (Images(Parameters(group="Fixed geometry", source="Images/Magnetic/FluxTube/CuboidParallelFlux.png")),
          Coordsys(
            extent=[-100, -100; 100, 100],
            grid=[2, 2],
            component=[20, 20]),
          Window(
            x=0.16,
            y=0.15,
            width=0.6,
            height=0.6),
          Icon,
          Diagram,
        Documentation(info="<HTML>
<p>
Please refer to the description of  the enclosing sub-package <a href=\"FixedShape\">FixedShape</a> for a description of all elements of this package.
</p>
</HTML>"));
      equation
        A = a * b;
        G_m = (my_0 * my_r * A)/ l;

      end Cuboid;

    end FixedShape;

    extends Modelica.Icons.Library;

  annotation (Documentation(info="<HTML>
<p>
Please have a look at <a href=\"UsersGuide.FluxTubeConcept\">Flux tube concept</a> and <a href=\"UsersGuide.ReluctanceForceCalculation\">Reluctance forces</a> in the Users Guide for comments on the underlying concept, on the calculation of reluctance forces from lumped magnetic network models and on the thereof resulting flux tube categories and sub-packages.
</p>
<p>
Although either the magnetic reluctance <i>R<sub>m</sub></i> or its reciprocal, the magnetic permeance <i>G<sub>m</sub></i> both represent geometric and material properties of a flux tube, both quantities are computed for each flux tube of this package. This enables for convenient calculation of the net reluctance or net permeance of a magnetic circuit, e.g. for calculation of a device's inductance. In case of a parallel connection of flux tube elements, the individual permeances sum up to the net permeance. For a series connection of flux tube elements, the individual reluctances sum up to the net reluctance.
</p>
</HTML>"));

    package Force
    "Flux tubes with reluctance force generation; constant permeability"
      partial model PartialForce
      "Base class for flux tubes with reluctance force generation; constant permeability"

        extends Modelica_Magnetic.Interfaces.TwoPortComponent;

        parameter SI.RelativePermeability my_r = 1
        "Relative magnetic permeability";

        SI.Force F_m "Reluctance force";

        SI.Reluctance R_m "Magnetic reluctance";
        SI.Permeance G_m "Magnetic permeance";
        Real dGmBydx
        "Derivative of permeance with respect to armature position";
        parameter Integer dlBydx = 1
        "Derivative of flux tube's varying dimension with respect to armature position; set to +1 or -1";

        annotation (
          Coordsys(
            extent=[-100, -100; 100, 100],
            grid=[2, 2],
            component=[20, 20]),
          Window(
            x=0.16,
            y=0.15,
            width=0.6,
            height=0.6),
          Icon(
            Rectangle(extent=[-70,30; 70,-30], style(
                color=45,
                rgbcolor={255,128,0},
                fillColor=7,
                rgbfillColor={255,255,255},
                fillPattern=1)),
            Line(points=[-70,0; -90,0], style(
                color=45,
                rgbcolor={255,128,0},
                fillColor=7,
                rgbfillColor={255,255,255},
                fillPattern=1)),
            Line(points=[70,0; 90,0], style(
                color=45,
                rgbcolor={255,128,0},
                fillColor=7,
                rgbfillColor={255,255,255},
                fillPattern=1))),
          Diagram,
        Documentation(info="<HTML>
<p>
Please refer to the description of  the enclosing sub-package <a href=\"Force\">Force</a> for a description of all elements of this package.
</p>
</HTML>"));

      Modelica.Mechanics.Translational.Interfaces.Flange_b flange
        "Generated reluctance force at armature position"
        annotation (extent=[-10,-10; 10,10]);
      equation

        V_mag = Phi * R_m;
        flange.f = -F_m;

        R_m = 1/G_m;

        F_m = 0.5 * V_mag^2 * dGmBydx;

      end PartialForce;

      model HollowCylinderAxialFlux
      "(Hollow) cylinder with axial flux; constant permeability"

        extends Modelica_Magnetic.FluxTube.Force.PartialForce;

        SI.Length l = flange.s "Axial length (in direction of flux)"
                                                  annotation(Dialog(group="Variable geometry"));
        parameter SI.Radius r_i = 0 "Inner radius of (hollow) cylinder";
        parameter SI.Radius r_o = 0.01 "Outer radius of (hollow) cylinder";

        annotation (Images(Parameters(group="Variable geometry", source="Images/Magnetic/FluxTube/HollowCylinderAxialFlux.png")),
            Documentation(info="<html>
<p>
Please refer to the description of  the enclosing sub-package <a href=\"Force\">Force</a> for a description of all elements of this package.
</p>
</html>"));

        SI.MagneticFluxDensity B "Homogeneous flux density";

    protected
        parameter SI.Area A = pi*(r_o^2 - r_i^2)
        "Cross-sectional area orthogonal to direction of flux";

      equation
        G_m = my_0*my_r * A /l;

        dGmBydx = -1 * my_0*my_r * A /l^2 * dlBydx;

        B = Phi/A;

      end HollowCylinderAxialFlux;
      annotation (Documentation(info="<html>
<p>
Please have a look at <a href=\"UsersGuide.ReluctanceForceCalculation\">Reluctance forces</a> in the Users Guide for an explanation of the different flux tube categories and resulting sub-packages.
</p>
Flux tube elements with generation of a reluctance force are intended for modelling of position-dependent air gap sections and permanent magnet sections respectively of translatory actuators. By default, the position co-ordinate of the mechanical connector flange.s is identical with the dimension l of the package's flux tube elements. l is the dimension changes with armature motion. If needed, the identity l=flange.s can be replaced by an actuator-specific equation, for example, when a flux tube length increases with decreasing armature position. The position co-ordinate of an element's translatory connector flange.s in turn will be identical with the armature position x in most cases, as the examples illustrate. </p>
<p>
The derivative of each element's permeance with respect to armature position dGmBydx is calculated from the derivative of the flux tube's permeance with respect to its varying dimension dGmBydl and the derivative of this dimension with respect to armature position dlBydx:</p>
<pre>
    dGm   dGm   dl
    --- = --- * --
     dx    dl   dx
</pre>
</p>
<p>
The parameter dlBydx must be set in each flux tube element to +1 or -1 according to the definition of the armature co-ordinate and the position of the element in a device's magnetic circuit. Proper match between armature motion and resulting variation of the flux tube length assures that the element's reluctance force acts in the right direction.
</p>
<p>
The shapes of the flux tubes defined in this package are rather simple. Only one dimenion varies with armature motion. Flux tubes with more complex variations of dimensions with armature motion can be defined by extending the base class <a href=\"PartialForce\">PartialForce</a>, if needed. Determination of the analytic derivative dGmBydl could become more complex for those flux tubes.
</p>
</html>"));

      model HollowCylinderRadialFlux
      "Hollow cylinder with radial flux; constant permeability"

        extends Modelica_Magnetic.FluxTube.Force.PartialForce;

        SI.Length l = flange.s "Axial length (orthogonal to direction of flux)"
                                                  annotation(Dialog(group="Variable geometry"));
        parameter SI.Radius r_i = 0.01 "Inner radius of hollow cylinder";
        parameter SI.Radius r_o = 0.015 "Outer radius of hollow cylinder";

        annotation (Images(Parameters(group="Variable geometry", source="Images/Magnetic/FluxTube/HollowCylinderRadialFlux_force.png")),
            Documentation(info="<html>
<p>
Please refer to the description of  the enclosing sub-package <a href=\"Force\">Force</a> for a description of all elements of this package.
</p>
</html>"));

        SI.MagneticFluxDensity B_avg
        "Average flux density (at arithmetic mean radius)";

    protected
        SI.Area A_avg
        "Average cross-sectional area orthogonal to direction of flux (at arithmetic mean radius)";

      equation
        G_m = my_0*my_r * 2 * pi * l /Modelica.Math.log(r_o/r_i);

        dGmBydx = my_0*my_r * 2 * pi/Modelica.Math.log(r_o/r_i) * dlBydx;

        A_avg = pi*(r_i + r_o) * l;
        B_avg = Phi/A_avg;

      end HollowCylinderRadialFlux;

      model CuboidParallelFlux
      "Cuboid with flux in direction of motion, e.g. air gap with rectangular cross-section; constant permeability"

        extends Modelica_Magnetic.FluxTube.Force.PartialForce;

        SI.Length l = flange.s "Axial length (in direction of flux)"
                                                  annotation(Dialog(group="Variable geometry"));
        parameter SI.Breadth a = 0.01 "Breadth of rectangular cross-section";
        parameter SI.Height b = 0.01 "Height of rectangular cross-section";

        annotation (Images(Parameters(group="Variable geometry", source="Images/Magnetic/FluxTube/CuboidParallelFlux.png")),
            Documentation(info="<html>
<p>
Please refer to the description of  the enclosing sub-package <a href=\"Force\">Force</a> for a description of all elements of this package.
</p>
</html>"));

        SI.MagneticFluxDensity B "Homogeneous flux density";

    protected
        parameter SI.Area A = a*b
        "Cross-sectional area orthogonal to direction of flux";

      equation
        G_m = my_0*my_r * A /l;

        dGmBydx = -1 * my_0*my_r * A /l^2 * dlBydx;

        B = Phi/A;

      end CuboidParallelFlux;

      model CuboidOrthogonalFlux
      "Cuboid with flux orthogonal to direction of motion; constant permeability"

        extends Modelica_Magnetic.FluxTube.Force.PartialForce;

        SI.Length l = flange.s
        "Axial length (in direction of motion, orthogonal to flux)"
                                                  annotation(Dialog(group="Variable geometry"));
        parameter SI.Breadth a = 0.01 "Breadth of rectangular cross-section";
        parameter SI.Height b = 0.01 "Height of rectangular cross-section";

        annotation (Images(Parameters(group="Variable geometry", source="Images/Magnetic/FluxTube/CuboidOrthogonalFlux.png")),
            Documentation(info="<html>
<p>
Please refer to the description of  the enclosing sub-package <a href=\"Force\">Force</a> for a description of all elements of this package.
</p>
</html>"));

        SI.MagneticFluxDensity B "Homogeneous flux density";

    protected
        SI.Area A "Cross-sectional area orthogonal to direction of flux";

      equation
        A = a*l;
        G_m = my_0*my_r * A /b;

        dGmBydx = my_0*my_r * a /b * dlBydx;

        B = Phi/A;

      end CuboidOrthogonalFlux;

      model LeakageAroundPoles
      "Leakage flux tube around cylindrical or prismatic poles"

        extends Modelica_Magnetic.FluxTube.Force.PartialForce;
        SI.Length l = flange.s "Axial length (in direction of flux)"
          annotation(Dialog(group="Variable geometry"));

        annotation (Images(Parameters(group="Variable geometry", source="Images/Magnetic/FluxTube/LeakageAroundPoles.png")),
            Documentation(info="<html>
<p>
Please refer to the description of  the enclosing sub-package <a href=\"Force\">Force</a> for a description of all elements of this package.
</p>
<p>
Leakage flux around a prismatic or cylindric air gap between to poles can be described with this model. Due to its constant radius of the leakage field r_leak, the model is rather simple. Whereas in reality the leakage radius is approximately constant for air gap lengths l greater than this radius, it decreases with air gap lengths less than the leakage radius. This decrease for small air gaps is neglected here, since the influence of the leakage flux tube compared to that of the enclosed main air gap (connected in parallel) decreases for decreasing air gap length l.
</p>
</html>"));

        parameter SI.Length t = 0.1
        "Depth orthogonal to flux; mean circumference of flux tube in case of cylindrical poles";
        parameter SI.Radius r_leak = 0.01 "Radius of leakage field";

      equation
        G_m = 2 * my_0 * t /pi * Modelica.Math.log(1 + pi/2 * r_leak/l);

      //derivative at full length:
      //  dGmBydx = 2 * my_0 * t /pi * 1/(1 + pi/2 * r_leak/l) * (-1)*pi/2*r_leak/l^2  * dlBydx;
      //simplified:
        dGmBydx = - my_0 * t * r_leak * dlBydx / (l^2 *(1 + pi/2 * r_leak/l));

      end LeakageAroundPoles;

    end Force;

    package Leakage
    "Leakage flux tubes with position-independent permeance and hence no force generation; my_r=1"
      partial model PartialLeakage
      "Base class for leakage flux tubes with position-independent permeance and hence no force generation; my_r=1"

        extends Modelica_Magnetic.Interfaces.TwoPortComponent;

        SI.Reluctance R_m "Magnetic reluctance";
        SI.Permeance G_m "Magnetic permeance";

        annotation (
          Coordsys(
            extent=[-100, -100; 100, 100],
            grid=[2, 2],
            component=[20, 20]),
          Window(
            x=0.16,
            y=0.15,
            width=0.6,
            height=0.6),
          Icon(
            Rectangle(extent=[-70,30; 70,-30], style(
                color=45,
                rgbcolor={255,128,0},
                fillColor=7,
                rgbfillColor={255,255,255},
                fillPattern=1)),
            Line(points=[-70,0; -90,0], style(
                color=45,
                rgbcolor={255,128,0},
                fillColor=7,
                rgbfillColor={255,255,255},
                fillPattern=1)),
            Line(points=[70,0; 90,0], style(
                color=45,
                rgbcolor={255,128,0},
                fillColor=7,
                rgbfillColor={255,255,255},
                fillPattern=1))),
          Diagram,
        Documentation(info="<HTML>
<p>
Please refer to the description of  the enclosing sub-package <a href=\"Leakage\">Leakage</a> for a description of all elements of this package.
</p>
</HTML>"));

      equation
        V_mag = Phi * R_m;
        R_m = 1/G_m;

      end PartialLeakage;

      annotation (Documentation(info="<html>
<p>
Please have a look at <a href=\"UsersGuide.ReluctanceForceCalculation\">Reluctance forces</a> in the Users Guide for an explanation of the different flux tube categories and resulting sub-packages.
</p>
<p>
Except for the element <a href=\"LeakageWithCoefficient\">LeakageWithCoefficient</a>, the permeances of all elements of this package are calculated from their geometry. These flux tube elements are intended for modelling of leakage fields through vacuum, air and other media with a relative permeability my_r=1.
</p>
<p>
All dimensions are defined as parameters. As a result, the shape of these elements will remain constant during dynamic simulation of actuators and reluctance forces will not be generated in these flux tube elements. A simple leakage flux tube with reluctance force generation is provided with the element <a href=\"Force.LeakageAroundPoles\">Force.LeakageAroundPoles</a>. In cases where the accuracy of that element is not sufficient, the leakage elements of this package can be adapted and extended so that they are able to change their shape with armature motion and to generate reluctance forces. This requires an extension of the partial model <a href=\"Force.PartialForce\">Force.PartialForce</a>, a higher variability of the variables representing the flux tube's dimensions, definition of a relationship between armature position and these dimensions and determination of the analytic derivative dG_m/dx of the flux tube's permeance G_m with respect to armature position x.
</p>
</html>"));

      model QuarterCylinder
      "Leakage flux from one edge to the opposite plane through a quarter cylinder"

        extends PartialLeakage;

        parameter SI.Thickness t = 0.1
        "Depth orthogonal to flux (=2*pi*r for cylindrical pole and r>>distance between edge and plane)";

        annotation (Images(Parameters(source="Images/Magnetic/FluxTube/Leakage/QuarterCylinder.png")),
            Documentation(info="<html>
<p>
Please refer to the description of  the enclosing sub-package <a href=\"Leakage\">Leakage</a> for a description of all elements of this package.
</p>
</html>"));

      equation
        G_m = my_0 * 0.52 * t;

      end QuarterCylinder;

      model QuarterHollowCylinder
      "Leakage flux in circumferential direction through a quarter hollow cylinder"

        extends PartialLeakage;

        parameter SI.Thickness t = 0.1
        "Depth orthogonal to flux (=2*pi*r for cylindrical pole and r>>r_i)";
        parameter Real ratio = 1 "Constant ratio b/r_i";

        annotation (Images(Parameters(source="Images/Magnetic/FluxTube/Leakage/QuarterHollowCylinder.png")),
            Documentation(info="<html>
<p>
Please refer to the description of  the enclosing sub-package <a href=\"Leakage\">Leakage</a> for a description of all elements of this package.
</p>
</html>"));

      equation
        G_m = 2* my_0 * t * Modelica.Math.log(1 + ratio) /pi;

      end QuarterHollowCylinder;

      model HalfCylinder "Leakage flux through the edges of a half cylinder"

        extends PartialLeakage;

        parameter SI.Thickness t = 0.1
        "Depth orthogonal to flux (=2*pi*r for cylindrical pole and r>>distance between edges)";

        annotation (Images(Parameters(source="Images/Magnetic/FluxTube/Leakage/HalfCylinder.png")),
            Documentation(info="<html>
<p>
Please refer to the description of  the enclosing sub-package <a href=\"Leakage\">Leakage</a> for a description of all elements of this package.
</p>
</html>"));

      equation
        G_m = my_0 * 0.26 * t;

      end HalfCylinder;

      model HalfHollowCylinder
      "Leakage flux in circumferential direction through a half hollow cylinder"

        extends PartialLeakage;

        parameter SI.Thickness t = 0.1
        "Depth orthogonal to flux (=2*pi*r for cylindrical pole and r>>r_i)";
        parameter Real ratio = 1 "Constant ratio b/r_i";

        annotation (Images(Parameters(source="Images/Magnetic/FluxTube/Leakage/HalfHollowCylinder.png")),
            Documentation(info="<html>
<p>
Please refer to the description of  the enclosing sub-package <a href=\"Leakage\">Leakage</a> for a description of all elements of this package.
</p>
</html>"));

      equation
        G_m = my_0 * t * Modelica.Math.log(1 + ratio) /pi;

      end HalfHollowCylinder;

      model QuarterSphere
      "Leakage flux through the corners of a quarter sphere"

        extends PartialLeakage;

        parameter SI.Diameter d = 0.01 "Diameter of quarter sphere";

        annotation (Images(Parameters(source="Images/Magnetic/FluxTube/Leakage/QuarterSphere.png")),
            Documentation(info="<html>
<p>
Please refer to the description of  the enclosing sub-package <a href=\"Leakage\">Leakage</a> for a description of all elements of this package.
</p>
</html>"));

      equation
        G_m = my_0 * 0.077 * d;

      end QuarterSphere;

      model QuarterHollowSphere
      "Leakage flux through the edges of a qarter hollow sphere"

        extends PartialLeakage;

        parameter SI.Thickness t = 0.01 "Thickness of sperical shell";

        annotation (Images(Parameters(source="Images/Magnetic/FluxTube/Leakage/QuarterHollowSphere.png")),
            Documentation(info="<html>
<p>
Please refer to the description of  the enclosing sub-package <a href=\"Leakage\">Leakage</a> for a description of all elements of this package.
</p>
</html>"));

      equation
        G_m = my_0 * 0.25 * t;

      end QuarterHollowSphere;

      model EighthOfSphere
      "Leakage flux through one edge and the opposite plane of an eighth of a sphere"

        extends PartialLeakage;

        parameter SI.Radius r = 0.01 "Radius of eighth of sphere";

        annotation (Images(Parameters(source="Images/Magnetic/FluxTube/Leakage/EighthOfSphere.png")),
            Documentation(info="<html>
<p>
Please refer to the description of  the enclosing sub-package <a href=\"Leakage\">Leakage</a> for a description of all elements of this package.
</p>
</html>"));

      equation
        G_m = my_0 * 0.308 * r;

      end EighthOfSphere;

      model EighthOfHollowSphere
      "Leakage flux through one edge and the opposite plane of an eighth of a hollow sphere"

        extends PartialLeakage;

        parameter SI.Thickness t = 0.01 "Thickness of sperical shell";

        annotation (Images(Parameters(source="Images/Magnetic/FluxTube/Leakage/EighthOfHollowSphere.png")),
            Documentation(info="<html>
<p>
Please refer to the description of  the enclosing sub-package <a href=\"Leakage\">Leakage</a> for a description of all elements of this package.
</p>
</html>"));

      equation
        G_m = my_0 * 0.5 * t;

      end EighthOfHollowSphere;

      model CoaxCylindersEndFaces
      "Leakage flux between the end planes of a inner solid cylinder and a coaxial outer hollow cylinder"

        extends PartialLeakage;

        parameter SI.Radius r_0 = 10e-3 "Radius of inner solid cylinder";
        parameter SI.Radius r_1 = 17e-3 "Inner radius of outer hollow cylinder";
        parameter SI.Radius r_2 = 20e-3 "Outer radius of outer hollow cylinder";

        annotation (Images(Parameters(source="Images/Magnetic/FluxTube/Leakage/CoaxCylindersEndFaces.png")),
            Documentation(info="<html>
<p>
Please refer to the description of  the enclosing sub-package <a href=\"Leakage\">Leakage</a> for a description of all elements of this package.
</p>
</html>"));

      equation
        assert(Modelica.Math.log(r_2/r_1) >= 2*(r_2-r_1)/(r_1+r_0), "No proper values assigned to Radii r_0...r_2!");
        G_m = 2*my_0 * sqrt( ((r_1+r_0)/2 * Modelica.Math.log(r_2/r_1))^2 - (r_2-r_1)^2);

      end CoaxCylindersEndFaces;

      model LeakageWithCoefficient
      "Leakage reluctance with respect to the reluctance of a useful flux path (not for dynamic simulation of actuators)"

        extends PartialLeakage;

        parameter SI.LeakageCoefficient c_leak = 0.3
        "Ratio leakage flux/(leakage flux + useful flux) = leakage flux/total flux";

        annotation (Images(Parameters(group="Reference reluctance", source="Images/Magnetic/FluxTube/Leakage/LeakageWithCoefficient.png")),
            Documentation(info="<html>
<p>
Please refer to the description of  the enclosing sub-package <a href=\"Leakage\">Leakage</a> for a description of all elements of this package.
</p>
<p>
Differently from the other flux tube elements of this package that are calculated from their geometry, this leakage reluctance is calculated with reference to the total reluctance of a useful flux path. Please refer to the <b>Parameters</b> section for an illustration of the resulting magnetic network. Exploiting <i>Kirchhoff</i>'s generalized current law, the leakage reluctance is calculated by means of a leakage coefficient c_leak.
</p>
<p>
<b>Attention:</b> <br>
This element must <b>not</b> be used <b>for dynamic simulation of</b> electro-magneto-mechanical <b>actuators</b>, where the shape of at least one flux tube element with reluctance force generation in the useful flux path changes with armature motion (e.g. air gap). This change results in a non-zero derivative dG_m/dx of those elements permeance G_m with respect to armature position x, which in turn will lead to a non-zero derivative of the leakage element's permeance with respect to armature position. This would result in a reluctance force generated by the leakage element that is not accounted for properly. Instead, use the leakage element for static analyses of magnetic networks only as illustrated in <a href=\"Examples.ElectrodynamicActuator.MagneticCircuitModel\">Examples.ElectrodynamicActuator.MagneticCircuitModel</a>.
</p>
</html>"));

        SI.Reluctance R_mUsefulTot = Modelica.Constants.inf
        "Total reluctance of useful flux path as reference"
           annotation(Dialog(group="Reference reluctance"));
        //Default value to ensure proper assignment of reference reluctance in magnetic network model

      equation
        assert(R_mUsefulTot < Modelica.Constants.inf, "No proper value assigned to R_mUsefulTot!");
        c_leak * R_m = R_mUsefulTot * (1 - c_leak);   // Generalized Kirchhoff's current law

      end LeakageWithCoefficient;
    end Leakage;

    model ConstantReluctance "Constant reluctance"

      extends Modelica_Magnetic.Interfaces.TwoPortComponent;

      parameter SI.Reluctance R_m = 1 "Magnetic reluctance";

      annotation (
        Coordsys(
          extent=[-100, -100; 100, 100],
          grid=[2, 2],
          component=[20, 20]),
        Window(
          x=0.16,
          y=0.15,
          width=0.6,
          height=0.6),
        Icon(
          Rectangle(extent=[-70,30; 70,-30], style(
              color=45,
              rgbcolor={255,128,0},
              fillColor=7,
              rgbfillColor={255,255,255},
              fillPattern=1)),
          Line(points=[-70,0; -90,0], style(
              color=45,
              rgbcolor={255,128,0},
              fillColor=7,
              rgbfillColor={255,255,255},
              fillPattern=1)),
          Line(points=[70,0; 90,0], style(
              color=45,
              rgbcolor={255,128,0},
              fillColor=7,
              rgbfillColor={255,255,255},
              fillPattern=1))),
        Diagram,
      Documentation(info="<HTML>
<p>
This constant reluctance is provided for test purposes and simple magnetic network models. The reluctance is not calculated from geometry and permeability of a flux tube, but is provided as a parameter.
</p>
</HTML>"));

    equation
      V_mag = Phi * R_m;

    end ConstantReluctance;
  end FluxTube;


  package Sources
  "Sources of different complexity of magnetomotive force and magnetic flux"
    extends Modelica.Icons.Library;

    annotation (Coordsys(
        extent=[-100, 100; 100, -100],
        grid=[2, 2],
        component=[20, 20]), Window(
        x=0.45,
        y=0.01,
        width=0.35,
        height=0.49,
        library=1,
        autolayout=1),
    Documentation(info="<html>
<p>
This package contains sources of magnetic potential difference:
</p>
</html>"));

    model ConstantMMF "Constant magnetomotive force"

      extends Modelica_Magnetic.Interfaces.TwoPortComponent;
      parameter SI.MagnetomotiveForce theta "Magnetomotive force";

      annotation (
        Coordsys(
          extent=[-100, -100; 100, 100],
          grid=[2, 2],
          component=[20, 20]),
        Icon(
          Ellipse(extent=[-40, -40; 40, 40], style(color=45)),
          Line(points=[100, 0; 40, 0], style(color=45)),
          Line(points=[-40, 0; -100, 0], style(color=45)),
          Line(points=[-12, 14; 12, 14], style(color=45)),
          Line(points=[0, 14; 0, -14], style(color=45)),
          Line(points=[-12, -14; 12, -14], style(color=45))),
        Diagram,
        Window(
          x=0.48,
          y=0.25,
          width=0.6,
          height=0.6),
        Documentation(info="<html>
<p>
Magnetic circuits under steady-state conditions, i.e. with stationary magnetic fields (change of magnetic flux  d&Phi;/dt=0) can be described with constant sources of a magnetomotive force (mmf). Constant mmf's are imposed by
<ul>
<li>
coils with stationary current (di/dt=0) and </li>
<li>
permanent magnets modelled with <i>Thévenin</i>'s equivalent magnetic circuit. </li>
</ul>
For modelling of reluctance actuators with constant mmf sources it is assumed that the armature is fixed so that no motion-induced flux change d&Phi;/dt can occur.
</p>
</html>"));
    equation
      V_mag = theta;
    end ConstantMMF;

    model SignalMMF "Signal-controlled magnetomotive force"

      extends Modelica_Magnetic.Interfaces.TwoPortComponent;
      Modelica.Blocks.Interfaces.RealInput theta "Magnetomotive force"
        annotation (extent=[-10, -50; 10, -30], rotation=90);

      annotation (
        Coordsys(
          extent=[-100, -100; 100, 100],
          grid=[2, 2],
          component=[20, 20]),
        Window(
          x=0.27,
          y=0.04,
          width=0.61,
          height=0.73),
        Icon(
          Ellipse(extent=[-30, 30; 30, -30], style(color=45)),
          Line(points=[-100, 0; -30, 0], style(color=45)),
          Line(points=[30, 0; 100, 0], style(color=45)),
          Line(points=[0, 10; 0, -10], style(color=45)),
          Line(points=[-10, -10; 10, -10], style(color=45)),
          Line(points=[-10, 10; 10, 10], style(color=45))),
      Diagram,
        Documentation(info="<html>
<p>
In electromagnetic devices, a change of a coil's magnetic flux linkage &Psi; reacts on the electrical subsystem in that a voltage v is induced due to <i>Faraday</i>'s law:
<pre>
    v = - d&Psi;/dt
</pre>
This reaction can possibly be neglected for
<ul>
<li>
modelling ofelectromagnetic actuators under quasi-stationary conditions (slow current change, slow armature motion),
<li>
modelling of current-controlled electromagnetic actuators (ideal current source) and</li>
<li>
for system simulation where the system dynamics is not governed by an electromagnetic actuator, but by the surrounding subsystems.
</ul>
In these cases, the magnetomotive force (mmf) imposed by a coil can easily be modelled with a signal-controlled mmf source. Except for the neglected dynamics, steady-state actuator forces will be calculated properly in actuator models based on these sources.
</p>
</html>"));
    equation
      V_mag = theta;
    end SignalMMF;

    model ElectroMagneticConverter "Electro-magnetic energy conversion"

      Interfaces.PositiveMagneticPort p_mag annotation (extent=[90, 50; 110, 70]);
      Interfaces.NegativeMagneticPort n_mag annotation (extent=[90, -70; 110, -50]);
      Modelica.Electrical.Analog.Interfaces.PositivePin p_el
        annotation (extent=[-110, 50; -90, 70]);
      Modelica.Electrical.Analog.Interfaces.NegativePin n_el
        annotation (extent=[-110, -70; -90, -50]);
      SI.Voltage v "Voltage";
      SI.Current i "Current";
      SI.MagneticPotentialDifference V_mag "Magnetomotive force";
      SI.MagneticFlux Phi "Magnetic flux coupled into magnetic circuit";

      Real w = 1 "Number of turns" annotation(Dialog(group="Variables"));
      SI.CouplingCoefficient c_coupl = 1
      "Ratio of coil's complete flux linkage to flux linked with the magnetic circuit; 0 < c_coupl <= 1; 1...complete flux linked with magnetic circuit"
                                                                                    annotation(Dialog(group="Variables"));

      //for information only:
      SI.MagneticFlux Psi "Flux linkage for information only";
      SI.Inductance L_stat
      "Static inductance abs(Psi/i) for information only (Caution: L(i=0) set to zero!)";

      annotation (
        Coordsys(
          extent=[-100, -100; 100, 100],
          grid=[2, 2],
          component=[20, 20]),
        Diagram(
          Polygon(points=[-134, 63; -124, 60; -134, 57; -134, 63], style(
              color=9,
              fillColor=9,
              fillPattern=1)),
          Line(points=[-150, 60; -125, 60], style(color=9, fillColor=9)),
          Polygon(points=[141, -57; 151, -60; 141, -63; 141, -57], style(
              color=9,
              fillColor=9,
              fillPattern=1)),
          Line(points=[125, -60; 150, -60], style(color=9, fillColor=9)),
          Text(
            extent=[128, -56; 144, -41],
            string="Phi",
            style(color=9, fillColor=9)),
          Text(
            extent=[128, 64; 145, 79],
            string="Phi",
            style(
              color=9,
              gradient=2,
              fillColor=9)),
          Line(points=[-150, -59; -125, -59], style(color=9, fillColor=9)),
          Polygon(points=[-140, -56; -150, -59; -140, -62; -140, -56], style(
              color=9,
              fillColor=9,
              fillPattern=1)),
          Text(
            extent=[-141, -56; -124, -41],
            string="i",
            style(color=9, fillColor=9)),
          Text(
            extent=[-150, 63; -133, 78],
            string="i",
            style(color=9, fillColor=9)),
          Line(points=[124, 61; 149, 61], style(color=9, fillColor=9)),
          Polygon(points=[134, 64; 124, 61; 134, 58; 134, 64], style(
              color=9,
              fillColor=9,
              fillPattern=1))),
        Window(
          x=0.26,
          y=0.14,
          width=0.58,
          height=0.58),
        Icon(
          Rectangle(extent=[-70,100; 70,-100], style(
              pattern=0,
              fillColor=7,
              rgbfillColor={255,255,255})),
          Ellipse(extent=[-50,0; -30,20]),
          Line(points=[-40,60; -40,40]),
          Ellipse(extent=[-50,20; -30,40]),
          Ellipse(extent=[-50,-20; -30,0]),
          Ellipse(extent=[-50,-40; -30,-20]),
          Line(points=[-40,-40; -40,-60]),
          Rectangle(extent=[-54,40; -40,-40],   style(color=7, fillColor=7)),
          Line(points=[-40,60; -92,60]),
          Line(points=[-40,-60; -90,-60]),
          Line(points=[0,100; -70,100],    style(pattern=2)),
          Line(points=[-70,100; -70,-100],     style(color=77, pattern=2)),
          Line(points=[0,-100; -70,-100],    style(pattern=2)),
          Line(points=[70,100; 0,100],    style(
              color=45,
              pattern=2,
              fillColor=45)),
          Line(points=[70,-100; 0,-100],    style(color=45, pattern=2)),
          Line(points=[70,100; 70,-100],     style(
              color=45,
              pattern=2,
              fillColor=45)),
          Ellipse(extent=[-4,-34; 64,34],   style(color=45)),
          Line(points=[30,-60; 30,-34],   style(color=45)),
          Line(points=[18,0; 42,0],   style(color=45)),
          Line(points=[42,10; 42,-12],   style(color=45)),
          Line(points=[30,34; 30,60],   style(color=45)),
          Line(points=[30,60; 100,60],   style(color=45)),
          Line(points=[30,-60; 90,-60],    style(color=45)),
          Text(
            extent=[-120,140; 120,100],
            string="%name",
            style(color=3, rgbcolor={0,0,255})),
          Line(points=[18,10; 18,-12],   style(color=45))),
      Documentation(info="<html>
<p>
The electro-magnetic energy conversion is given by <i>Ampere</i>'s law and <i>Faraday</i>'s law respectively:
<pre>
    V_mag = c_coupl * i*w
    w * der(&Phi;) =  -c_coupl * v
</pre>
V_mag is the magnetomotive force that is supplied to the connected magnetic circuit, &Phi; is the magnetic flux through the associated branch of this magnetic circuit. The negative sign of the induced voltage v is due to <i>Lenz</i>'s law. The coupling coefficient c_coupl denotes leakage: Only a portion of the coil's flux linkage &Psi; contributes to the useful magnetic flux &Phi; in the magnetic circuit. The influence of c_coupl can be interpreted as a leakage inductance connected in series with the inductance that is effective in the magnetic circuit. <br>
<br>
The flux linkage &Psi; and the static inductance L_stat = &Psi;/i are calculated for information only. Note that L_stat is set to zero in case of zero current. <br>
<br>
The variability of the number of turns w and the coupling coefficient c_coupl is not restricted to <i>parameter</i>, although both variables can be defined as parameters in many modelling cases. The extended variability allows for calculation of w from other variables (i.e. from non-parameters) in a model or for a variable coupling coefficient, e.g. a position-dependent one in an actuator model.
</p>
</html>"));

    equation
      v = p_el.v - n_el.v;
      0 = p_el.i + n_el.i;
      i = p_el.i;

      V_mag = p_mag.V_mag - n_mag.V_mag;
      0 = p_mag.Phi + n_mag.Phi;
      Phi = p_mag.Phi;

      assert(c_coupl>0 and c_coupl<=1, "c_coupl out of allowed range 0 <= c_coupl <= 1");

      //converter equations:
      V_mag = c_coupl * i*w;   // Ampere's law
      w*der(Phi) = - c_coupl * v;   // Faraday's law

      //for information only:
      Psi = w * Phi / c_coupl;
      //use of abs() for positive results; due to Modelica sign conventions for flow into connectors:
      L_stat = if i>0 or i<0 then abs(Psi/i) else 0;

    end ElectroMagneticConverter;

    model ConstantFlux "Source of constant magnetic flux"

      extends Modelica_Magnetic.Interfaces.TwoPortComponent;
      parameter SI.MagneticFlux Phi_source = 1 "Magnetic source flux";

      annotation (
        Coordsys(
          extent=[-100, -100; 100, 100],
          grid=[2, 2],
          component=[20, 20]),
        Icon(
          Ellipse(extent=[-40, -40; 40, 40], style(color=45)),
          Line(points=[100, 0; 40, 0], style(color=45)),
          Line(points=[-40, 0; -100, 0], style(color=45)),
          Line(points=[0,40; 0,-40], style(color=45, rgbcolor={255,128,0}))),
        Diagram(
          Line(points=[-125, 0; -115, 0], style(color=9)),
          Line(points=[-120, -5; -120, 5], style(color=9)),
          Line(points=[115, 0; 125, 0], style(color=9)),
          Ellipse(extent=[-40,-40; 40,40],   style(color=45)),
          Line(points=[90,0; 40,0],    style(color=45)),
          Line(points=[-40,0; -90,0],    style(color=45)),
          Line(points=[0,40; 0,-40], style(color=45, rgbcolor={255,128,0}))),
        Window(
          x=0.48,
          y=0.25,
          width=0.6,
          height=0.6),
        Documentation(info="<html>
<p>
Sources of a constant magnetic flux are useful for modelling of permanent magnets with <i>Norton</i>'s magnetic equivalent circuit.
</p>
</html>"));
    equation
      Phi = Phi_source;
    end ConstantFlux;

    model SignalFlux "Signal-controlled magnetic flux source"

      extends Modelica_Magnetic.Interfaces.TwoPortComponent;
      Modelica.Blocks.Interfaces.RealInput Phi_source "Imposed magnetic flux"
        annotation (extent=[-10, -50; 10, -30], rotation=90);

      annotation (
        Coordsys(
          extent=[-100, -100; 100, 100],
          grid=[2, 2],
          component=[20, 20]),
        Window(
          x=0.27,
          y=0.04,
          width=0.61,
          height=0.73),
        Icon(
          Ellipse(extent=[-30, 30; 30, -30], style(color=45)),
          Line(points=[-100, 0; -30, 0], style(color=45)),
          Line(points=[30, 0; 100, 0], style(color=45)),
          Line(points=[0,30; 0,-30], style(color=45, rgbcolor={255,128,0}))),
      Diagram,
        Documentation(info="<html>
<p>
This source of a magnetic flux is intended for test purposes, e.g. for simulation and subsequent plotting of a softmagnetic material's magnetisation characteristics if used together with a nonlinear reluctance element.
</p>
</html>"));
    equation
      Phi = Phi_source;
    end SignalFlux;

public
    model CoilDesign
    "Calculation of winding parameters (wire diameter, number of turns et al.) and recalculation with optionally chosen parameters; to be adapted to particular design tasks"

      parameter SI.Resistivity rho_20 = 0.0178e-6
      "Resistivity of conductor material at 20°C";     //default material: Copper
      parameter Utilities.TemperatureCoefficient alpha = 0.0039
      "Temperature coefficient of conductor material's resistivity";     //default material: Copper
      parameter SI.CelsiusTemperature T_opCelsius = 20
      "Winding's operating temperature";

      final parameter SI.Resistivity rho = rho_20 * (1 + alpha *(T_opCelsius - 20))
      "Resistivity at operating temperature";

      parameter SI.Height h_w "Height of winding's cross-section";
      parameter SI.Breadth b_w "Breadth of winding's cross-section";

      final parameter SI.Area A_w = h_w * b_w "Cross-section area of winding";

      parameter SI.Length l_avg "Average length of one turn";

      parameter SI.Voltage U
      "Operating voltage (nominal/ minimum/ maximum voltage depending on design objective)";

      parameter SI.CurrentDensity J_desired = 4e6
      "DESIRED current density at operating temperature and voltage resp.";

      parameter Real c_condFillChosen = 0.6
      "CHOSEN conductor filling factor = total conductor area without insulation/ total winding area";

      final parameter Real w_calculated = U/ (rho * l_avg * J_desired)
      "CALCULATED number of turns";

      final parameter SI.Diameter d_wireCalculated = sqrt(4 * A_w * c_condFillChosen /(pi * w_calculated))
      "CALCULATED wire diameter (without insulation)";

      final parameter SI.Area A_wireCalculated = pi * d_wireCalculated^2 / 4
      "Calculated wire cross-section area";

      final parameter SI.Resistance R_calculated = rho * w_calculated * l_avg / A_wireCalculated
      "Winding resistance at operating temperature and voltage resp. with CALCULATED number of turns and wire diameter";

      final parameter SI.Power P_calculated = U^2 / R_calculated
      "Winding's ohmic losses at operating temperature and voltage resp. with CALCULATED number of turns and wire diameter";

      parameter SI.Diameter d_wireChosen = d_wireCalculated
      "CHOSEN available wire diameter (without insulation)"   annotation(Dialog(group="Chosen feasible parameters (optional)"));

      parameter Real w_chosen = w_calculated "CHOSEN number of turns"
                                                       annotation(Dialog(group="Chosen feasible parameters (optional)"));

      final parameter SI.Area A_wireChosen = pi * d_wireChosen^2 / 4
      "Wire cross-section area resulting from CHOSEN wire diameter";

      final parameter SI.Resistance R_actual = rho * w_chosen * l_avg / A_wireChosen
      "Winding resistance at operating temperature and voltage resp. resulting from CHOSEN number of turns and wire diameter";

      final parameter SI.Power P_actual = U^2 / R_actual
      "Winding's ohmic losses at operating temperature and voltage resp. resulting from CHOSEN number of turns and wire diameter";

      final parameter SI.CurrentDensity J_actual = U * 4/(R_actual * pi * d_wireChosen^2)
      "Current density at operating temperature and voltage resp. resulting from CHOSEN number of turns and wire diameter";

      final parameter Real c_condFillActual = w_chosen * pi * d_wireChosen^2 /(4 * A_w)
      "Conductor filling factor resulting from CHOSEN number of turns and wire diameter";

      annotation (Icon(
          Rectangle(extent=[-100,100; 100,-100], style(
              color=45,
              rgbcolor={255,128,0},
              fillColor=7,
              rgbfillColor={255,255,255})),
          Text(
            extent=[-136,144; 136,104],
            style(color=3, rgbcolor={0,0,255}),
            string="%name"),
          Text(
            extent=[-84,80; 56,0],
            style(color=45, rgbcolor={255,128,0}),
            string="R(T)"),
          Text(
            extent=[-84,0; -24,-78],
            style(color=45, rgbcolor={255,128,0}),
            string="w"),
          Line(points=[-100,60; 100,60], style(color=45, rgbcolor={255,128,0})),
          Line(points=[-100,20; 100,20], style(color=45, rgbcolor={255,128,0})),
          Line(points=[-100,-20; 100,-20], style(color=45, rgbcolor={255,128,0})),
          Line(points=[-100,-60; 100,-60], style(color=45, rgbcolor={255,128,0})),
          Line(points=[-60,-100; -60,100], style(color=45, rgbcolor={255,128,0})),
          Line(points=[-20,-100; -20,100], style(color=45, rgbcolor={255,128,0})),
          Line(points=[20,-100; 20,100], style(color=45, rgbcolor={255,128,0})),
          Line(points=[60,-100; 60,100], style(color=45, rgbcolor={255,128,0}))),
                              Diagram,
        Documentation(info="<html>
<p>
This model exemplarily shows dimensioning of a winding (wire diameter, number of turns) based on desired operating conditions (voltage, temperature, current density, conductor filling factor) for a given cross-section area of the winding. It can be modified according to the parameters given and sought after for a particular design project. <br>
<br>
The calculated winding resistance and number of turns can be used as input parameters to the electrical subsystem
of a device to be modelled. Operating voltage U can be minimum, nominal and maximum voltage respectively as specified for a particular design project. In conjunction with the setting of the operating temperature T_opCelsius, this enables for analysis of the device under worst-case conditions (e.g. minimum required magnetomotive force, maximum allowed ohmic losses, minimum and maximum force respectively).<br>
<br>
For manufacturing of a winding, the obtained wire diameter d_wireCalculated must be rounded to that of an available wire. In order to analyse the influence of this rounding, one can enter the chosen wire diameter d_wireChosen and number of turns w_chosen as optional input. Calculation of the resulting winding parameters enables for comparison with the ones obtained otherwise. <br>
<br>
</p>
</html>"));

    end CoilDesign;

  end Sources;


  package Material
  "Magnetisation characteristics of common soft magnetic and hard magnetic materials"
    package SoftMagnetic
    "Characteristics my_r(B) of common soft magnetic materials; hysteresis neglected"
      package Steel "Various ferromagnetic steels"
        record Steel_9SMnPb28 =
          Modelica_Magnetic.Material.SoftMagnetic.ApproximationData (
            my_i=400,
            B_myMax=1.488,
            c_a=1200,
            c_b=3,
            n=12.5) "9SMnPb28 (1.0718)"   annotation (Documentation(info="<html>
<p>
Please refer to the description of  the enclosing package <a href=\"SoftMagnetic\">SoftMagnetic</a> for a description of all soft magnetic material characteristics of this package.
</p>
</html>"));
        record Steel_9SMn28K =
          Modelica_Magnetic.Material.SoftMagnetic.ApproximationData (
            my_i=500,
            B_myMax=1.036,
            c_a=43414,
            c_b=35.8,
            n=14) "9SMn28k (1.0715)"   annotation (Documentation(info="<html>
<p>
Please refer to the description of  the enclosing package <a href=\"SoftMagnetic\">SoftMagnetic</a> for a description of all soft magnetic material characteristics of this package.
</p>
</html>"));
        record DC01 =
          Modelica_Magnetic.Material.SoftMagnetic.ApproximationData (
            my_i=5,
            B_myMax=1.1,
            c_a=6450,
            c_b=3.65,
            n=7.7) "DC01 (1.0330, previously St2)"   annotation (
            Documentation(info="<html>
<p>
Please refer to the description of  the enclosing package <a href=\"SoftMagnetic\">SoftMagnetic</a> for a description of all soft magnetic material characteristics of this package.
</p>
</html>"));
        record DC03 =
          Modelica_Magnetic.Material.SoftMagnetic.ApproximationData (
            my_i=0,
            B_myMax=1.05,
            c_a=27790,
            c_b=16,
            n=10.4) "DC03 (1.0347, previously St3)"   annotation (
            Documentation(info="<html>
<p>
Please refer to the description of  the enclosing package <a href=\"SoftMagnetic\">SoftMagnetic</a> for a description of all soft magnetic material characteristics of this package.
</p>
</html>"));
        record X6Cr17 =
          Modelica_Magnetic.Material.SoftMagnetic.ApproximationData (
            my_i=274,
            B_myMax=1.1,
            c_a=970,
            c_b=1.2,
            n=8.3) "X6Cr17 (1.4016)"   annotation (Documentation(info="<html>
<p>
Please refer to the description of  the enclosing package <a href=\"SoftMagnetic\">SoftMagnetic</a> for a description of all soft magnetic material characteristics of this package.
</p>
</html>"));
        record AISI_1008 =
          Modelica_Magnetic.Material.SoftMagnetic.ApproximationData (
            my_i=200,
            B_myMax=1.17,
            c_a=8100,
            c_b=2.59,
            n=10) "AISI 1008 (1.0204)"   annotation (Documentation(info="<html>
<p>
Please refer to the description of  the enclosing package <a href=\"SoftMagnetic\">SoftMagnetic</a> for a description of all soft magnetic material characteristics of this package.
</p>
</html>"));
        record AISI_12L14 =
          Modelica_Magnetic.Material.SoftMagnetic.ApproximationData (
            my_i=10,
            B_myMax=0.94,
            c_a=5900,
            c_b=4.19,
            n=6.4) "AISI 12L14 (1.0718)"   annotation (Documentation(info="<html>
<p>
Please refer to the description of  the enclosing package <a href=\"SoftMagnetic\">SoftMagnetic</a> for a description of all soft magnetic material characteristics of this package.
</p>
</html>"));
        annotation (Documentation(info="<html>
<p>
Please refer to the description of  the enclosing package <a href=\"SoftMagnetic\">SoftMagnetic</a> for a description of all soft magnetic material characteristics of this package.
</p>
</html>"));
      end Steel;

      package ElectricSheet "Various electric sheets"
        extends Modelica.Icons.Library;

      record M330_50A =
        Modelica_Magnetic.Material.SoftMagnetic.ApproximationData (
            my_i=500,
            B_myMax=0.7,
            c_a=24000,
            c_b=9.38,
            n=9.6) "M330-50A (1.0809) @ 50Hz" annotation (Documentation(info="<html>
<p>
Please refer to the description of  the enclosing package <a href=\"SoftMagnetic\">SoftMagnetic</a> for a description of all soft magnetic material characteristics of this package.
</p>
<p>
Sample: complete core after machining and packet assembling<br>
</p>
</html>"));

      record M350_50A =
        Modelica_Magnetic.Material.SoftMagnetic.ApproximationData (
            my_i=1210,
            B_myMax=1.16,
            c_a=24630,
            c_b=2.44,
            n=14) "M350-50A (1.0810) @ 50Hz" annotation (Documentation(info="<html>
<p>
Please refer to the description of  the enclosing package <a href=\"SoftMagnetic\">SoftMagnetic</a> for a description of all soft magnetic material characteristics of this package.
</p>
<p>
Sample: sheet strip<br>
Measurement: Epstein frame
</p>
</html>"));

      record M530_50A =
        Modelica_Magnetic.Material.SoftMagnetic.ApproximationData (
            my_i=2120,
            B_myMax=1.25,
            c_a=12400,
            c_b=1.6,
            n=13.5) "M530-50A (1.0813) @ 50Hz" annotation (Documentation(
              info="<html>
<p>
Please refer to the description of  the enclosing package <a href=\"SoftMagnetic\">SoftMagnetic</a> for a description of all soft magnetic material characteristics of this package.
</p>
<p>
Sample: sheet strip<br>
Measurement: Epstein frame
</p>
</html>"));

      record M700_100A =
        Modelica_Magnetic.Material.SoftMagnetic.ApproximationData (
            my_i=1120,
            B_myMax=1.2,
            c_a=20750,
            c_b=3.55,
            n=13.15) "M700-100A (1.0826) @ 50Hz" annotation (Documentation(
              info="<html>
<p>
Please refer to the description of  the enclosing package <a href=\"SoftMagnetic\">SoftMagnetic</a> for a description of all soft magnetic material characteristics of this package.
</p>
<p>
Sample: sheet strip<br>
Measurement: Epstein frame
</p>
</html>"));

      record M940_100A =
        Modelica_Magnetic.Material.SoftMagnetic.ApproximationData (
            my_i=680,
            B_myMax=1.26,
            c_a=17760,
            c_b=3.13,
            n=13.9) "M940-100A @ 50Hz" annotation (Documentation(info="<html>
<p>
Please refer to the description of  the enclosing package <a href=\"SoftMagnetic\">SoftMagnetic</a> for a description of all soft magnetic material characteristics of this package.
</p>
<p>
Sample: sheet strip<br>
Measurement: Epstein frame
</p>
</html>"));

        annotation (Documentation(info="<html>
<p>
Please refer to the description of  the enclosing package <a href=\"SoftMagnetic\">SoftMagnetic</a> for a description of all soft magnetic material characteristics of this package.
</p>
</html>"));
      end ElectricSheet;

      package PureIron "Pure iron "
        record RFe80 =
          Modelica_Magnetic.Material.SoftMagnetic.ApproximationData (
            my_i=123,
            B_myMax=1.27,
            c_a=44410,
            c_b=6.4,
            n=10) "Hyperm 0 (RFe80)"    annotation (Documentation(info="<html>
<p>
Please refer to the description of  the enclosing package <a href=\"SoftMagnetic\">SoftMagnetic</a> for a description of all soft magnetic material characteristics of this package.
</p>
<p>
Source of B(H) characteristics: Product catalogue <i>Magnequench</i>, 2000
</p>
</html>"));
        record VacoferS2 =
          Modelica_Magnetic.Material.SoftMagnetic.ApproximationData (
            my_i=2666,
            B_myMax=1.15,
            c_a=187000,
            c_b=4.24,
            n=19) "VACOFER S2 (99.95% Fe)"    annotation (Documentation(info="<html>
<p>
Please refer to the description of  the enclosing package <a href=\"SoftMagnetic\">SoftMagnetic</a> for a description of all soft magnetic material characteristics of this package.
</p>
<p>
Source of B(H) characteristics:
<dd>
<p><i>Boll, R.</i>: Weichmagnetische Werkstoffe: Einf&uuml;hrung in den Magnetismus, VAC-Werkstoffe und ihre Anwendungen. 4th ed. Berlin, M&uuml;nchen: Siemens Aktiengesellschaft 1990</p>
</dd>
</p>
</html>"));
      end PureIron;

      extends Modelica.Icons.Library;

      package CobaltIron "Cobalt iron"
        record Vacoflux50 =
          Modelica_Magnetic.Material.SoftMagnetic.ApproximationData (
            my_i=3850,
            B_myMax=1.75,
            c_a=11790,
            c_b=2.63,
            n=15.02) "Vacoflux 50 (50% CoFe)"
                                           annotation (Documentation(info="<html>
<p>
Please refer to the description of  the enclosing package <a href=\"SoftMagnetic\">SoftMagnetic</a> for a description of all soft magnetic material characteristics of this package.
</p>
<p>
Source of B(H) characteristics: VACUUMSCHMELZE GmbH & Co. KG, Germany
</p>
</html>"));
      end CobaltIron;

      package NickelIron "Nickel iron"
        record MuMetall =
          Modelica_Magnetic.Material.SoftMagnetic.ApproximationData (
            my_i=27300,
            B_myMax=0.46,
            c_a=1037500,
            c_b=3.67,
            n=10) "MUMETALL (77% NiFe)"   annotation (Documentation(info="<html>
<p>
Please refer to the description of  the enclosing package <a href=\"SoftMagnetic\">SoftMagnetic</a> for a description of all soft magnetic material characteristics of this package.
</p>
<p>
Source of B(H) characteristics:
<dd>
<p><i>Boll, R.</i>: Weichmagnetische Werkstoffe: Einf&uuml;hrung in den Magnetismus, VAC-Werkstoffe und ihre Anwendungen. 4th ed. Berlin, M&uuml;nchen: Siemens Aktiengesellschaft 1990</p>
</dd>
</p>
</html>"));
        record Permenorm3601K3 =
          Modelica_Magnetic.Material.SoftMagnetic.ApproximationData (
            my_i=3000,
            B_myMax=0.67,
            c_a=50000,
            c_b=2.39,
            n=9.3) "PERMENORM 3601 K3 (36% NiFe)"    annotation (Documentation(
              info="<html>
<p>
Please refer to the description of  the enclosing package <a href=\"SoftMagnetic\">SoftMagnetic</a> for a description of all soft magnetic material characteristics of this package.
</p>
<p>
Source of B(H) characteristics:
<dd>
<p><i>Boll, R.</i>: Weichmagnetische Werkstoffe: Einf&uuml;hrung in den Magnetismus, VAC-Werkstoffe und ihre Anwendungen. 4th ed. Berlin, M&uuml;nchen: Siemens Aktiengesellschaft 1990</p>
</dd>
</p>
</html>"));
      end NickelIron;

      function my_rApprox
      "Approximation of relative permeability my_r as a function of flux density B for soft magnetic materials"

        extends Modelica.Icons.Function;

        input Real B "Flux density in ferromagnetic flux tube element";
      //Material specific parameter set:
        input Real my_i "Initial relative permeability at B=0";
        input Real B_myMax "Flux density at maximum relative permeability";
        input Real c_a "Coefficient of approximation function";
        input Real c_b "Coefficient of approximation function";
        input Real n "Exponent of approximation function";

        output Real my_r
        "Relative magnetic permeability of ferromagnetic flux tube element";

    protected
        Real B_N
        "Flux density B normalized to flux density at maximum relative permeability B_myMax";

      algorithm
        B_N := abs(B/B_myMax);
        my_r := 1 + (my_i-1 + c_a*B_N)/(1 + c_b*B_N + B_N^n);

        annotation (Documentation(info="<html>
<p>
The relative permeability my_r as a function of flux density B for all soft magnetic materials currently included in this library is approximated with the following function <a href=\"UsersGuide.Literature\">[4]</a>:
<br>
<dl>
<dd>
<p><IMG SRC=\"../Images/Magnetic/Material/SoftMagnetic/eq_my_rApprox.png\" ALT=\"Equation for approximation my_r(B)\"></p>
</dd>
</dl>
<br>
Two of the five parameters of this equation have a physical meaning, namely the initial relative permeability my_i at B=0 and the magnetic flux density at maximum permeability B_myMax. B_N is the flux density normalized to latter parameter.
</html>"));
      end my_rApprox;
      annotation (Documentation(info="<html>
<p>
The magnetisation characteristics my_r(B) of all soft magnetic materials currently included in this library are approximated with a <a href=\"my_rApprox\">function</a>. Each material is characterised by the five parameters of this function. The approximated characteristics my_r(B) for most of the ferromagnetic materials currently included are shown in the plots below (solid lines) together with the original data points compiled from measurements and literature.
<dl>
<dd>
<p><IMG SRC=\"../Images/Magnetic/Material/SoftMagnetic/Steel.png\" ALT=\"Approximated magnetization characteristics of selected steels\"></p>
</dd>
<dd>
<p><IMG SRC=\"../Images/Magnetic/Material/SoftMagnetic/Miscellaneous.png\" ALT=\"Approximated magnetization characteristics of miscellaneous soft magnetic materials\"></p>
</dd>
<dd>
<p><IMG SRC=\"../Images/Magnetic/Material/SoftMagnetic/ElectricSheet.png\" ALT=\"Approximated magnetization characteristics of included electric sheets\"></p>
</dd>
</dl>
</p>
<p>
For the nonlinear curve fit, data points for high flux densities (approximately B>1T) have been weighted higher than the ones for low flux densities. This is due to the large impact of saturated ferromagnetic sections in a magnetic circuit compared to that of non-saturated sections with relative permeabilities my_r>>1.
</p>
<p>
Note that the magnetisation characteristics largely depend on possible previous machining and on measurement conditions. A virgin material normally has a considerably higher permeability than the same material after machining (and packet assembling in case of electric sheets). This is indicated in the above plots by different magnetisation curves for similar materials. In most cases, the original data points represent commutating curves obtained with measurements at 50Hz.
</p>
<p>
Additional user-specific materials can be defined as needed. This requires determination of the approximation parameters from the original data points, preferably with a nonlinear curve fit.
</p>
</html>"));

      record ApproximationData
      "Coefficients for approximation of soft magnetic materials"

        extends Modelica.Icons.Record;

        parameter SI.RelativePermeability my_i = 1
        "Initial relative permeability at B=0";

        parameter SI.MagneticFluxDensity B_myMax = 1
        "Flux density at maximum relative permeability";

        parameter Real c_a = 1 "Coefficient of approximation function";

        parameter Real c_b = 1 "Coefficient of approximation function";

        parameter Real n = 1 "Exponent of approximation function";

        annotation (Documentation(info="<html>
<p>
The parameters needed for <a href=\"my_rApprox\">approximation of the magnetisation characteristics</a> of included soft magnetic materials are declared in this record.
</p>
</html>"));
      end ApproximationData;
    end SoftMagnetic;

    extends Modelica.Icons.Library;

  annotation (Documentation(info="<html>
</html>"));

    package HardMagnetic
    "Characteristics of common permanent magnetic materials (temperature dependence considered)"

      model PermanentMagnetBehaviour
      "Relative permeability and temperature-dependent coercivity of permanent magnetic materials"

        replaceable record material =
          Modelica_Magnetic.Material.HardMagnetic.PermanentMagnetData
        "Characteristics of common permanent magnetic materials"
                                                  annotation(choicesAllMatching=true, Dialog(group="Material"));

        material mat;

        parameter SI.CelsiusTemperature T_opCelsius = 20
        "Operating temperature";

        final parameter SI.MagneticFluxDensity B_r = mat.B_rRef * (1 + mat.alpha_Br *(T_opCelsius - mat.T_refCelsius))
        "Remanence at operating temperature";
        final parameter SI.MagneticFieldStrength H_cB = mat.H_cBRef * (1 + mat.alpha_Br *(T_opCelsius - mat.T_refCelsius))
        "Coercivity at operating temperature";
        final parameter SI.RelativePermeability my_r = B_r/ (my_0 * H_cB)
        "Relative permeability";

        annotation (Icon(
            Rectangle(extent=[-100,100; 100,-100], style(
                color=45,
                rgbcolor={255,128,0},
                fillColor=7,
                rgbfillColor={255,255,255})),
            Line(points=[-70,90; -60,100], style(
                color=45,
                rgbcolor={255,128,0},
                fillPattern=1)),
            Line(points=[-100,0; -60,40], style(
                color=45,
                rgbcolor={255,128,0},
                fillPattern=1)),
            Line(points=[-50,50; -40,60], style(
                color=45,
                rgbcolor={255,128,0},
                fillPattern=1)),
            Line(points=[-30,70; 0,100], style(
                color=45,
                rgbcolor={255,128,0},
                fillPattern=1)),
            Line(points=[-100,-60; -90,-50], style(
                color=45,
                rgbcolor={255,128,0},
                fillPattern=1)),
            Line(points=[-80,-40; -40,0], style(
                color=45,
                rgbcolor={255,128,0},
                fillPattern=1)),
            Line(points=[-30,10; -20,20], style(
                color=45,
                rgbcolor={255,128,0},
                fillPattern=1)),
            Line(points=[-10,30; 30,70], style(
                color=45,
                rgbcolor={255,128,0},
                fillPattern=1)),
            Line(points=[40,80; 50,90], style(
                color=45,
                rgbcolor={255,128,0},
                fillPattern=1)),
            Line(points=[-80,80; -100,60], style(
                color=45,
                rgbcolor={255,128,0},
                fillPattern=1)),
            Line(points=[-10,-30; 0,-20], style(
                color=45,
                rgbcolor={255,128,0},
                fillPattern=1)),
            Line(points=[10,-70; 20,-60], style(
                color=45,
                rgbcolor={255,128,0},
                fillPattern=1)),
            Line(points=[10,-10; 50,30], style(
                color=45,
                rgbcolor={255,128,0},
                fillPattern=1)),
            Line(points=[60,40; 70,50], style(
                color=45,
                rgbcolor={255,128,0},
                fillPattern=1)),
            Line(points=[100,80; 80,60], style(
                color=45,
                rgbcolor={255,128,0},
                fillPattern=1)),
            Line(points=[-60,-80; -20,-40], style(
                color=45,
                rgbcolor={255,128,0},
                fillPattern=1)),
            Line(points=[-80,-100; -70,-90], style(
                color=45,
                rgbcolor={255,128,0},
                fillPattern=1)),
            Line(points=[30,-50; 70,-10], style(
                color=45,
                rgbcolor={255,128,0},
                fillPattern=1)),
            Line(points=[80,0; 90,10], style(
                color=45,
                rgbcolor={255,128,0},
                fillPattern=1)),
            Line(points=[0,-80; -20,-100], style(
                color=45,
                rgbcolor={255,128,0},
                fillPattern=1)),
            Line(points=[50,-90; 90,-50], style(
                color=45,
                rgbcolor={255,128,0},
                fillPattern=1)),
            Text(
              extent=[-136,144; 136,104],
              style(color=3, rgbcolor={0,0,255}),
              string="%name"),
            Text(
              extent=[12,82; 98,2],
              style(color=45, rgbcolor={255,128,0}),
              string="(T)"),
            Text(
              extent=[-90,82; -30,2],
              style(color=45, rgbcolor={255,128,0}),
              string="H"),
            Text(
              extent=[-40,52; 10,-8],
              style(color=45, rgbcolor={255,128,0}),
              string="cB"),
            Text(
              extent=[-94,2; -34,-78],
              style(color=45, rgbcolor={255,128,0}),
              string="µ"),
            Text(
              extent=[-56,-28; -16,-88],
              style(color=45, rgbcolor={255,128,0}),
              string="r")),     Diagram,
          Documentation(info="<html>
<p>
Please refer to the description of  the enclosing package <a href=\"HardMagnetic\">HardMagnetic</a> for a description of all elements of this package.
</p>
<p>
In the records defining the characteristics of a permanent magnetic materials, remanence B_rRef and coercivity H_cBRef are given for a reference temperature T_refCelsius, usually 20°C. Using the also defined temperature coefficient of remanence alpha_Br, remanence B_r and coercivity H_cB are calculated for a given operating temperature of the permanent magnet T_opCelsius. In addition, the relative permeability my_r is calculated. Have a look at <a href=\"Examples.ElectrodynamicActuator.MagneticCircuitModel\">Examples.ElectrodynamicActuator.MagneticCircuitModel</a> for an exemplary use of this component.
</p>
</html>"));

      end PermanentMagnetBehaviour;

      record NdFeB =
        Modelica_Magnetic.Material.HardMagnetic.PermanentMagnetData (
          H_cBRef=900000,
          B_rRef=1.2,
          T_refCelsius=20,
          alpha_Br=-0.001) "NdFeB sintered; exemplary values"
                                      annotation (Documentation(info="<html>
<p>
Please refer to the description of  the enclosing package <a href=\"HardMagnetic\">HardMagnetic</a> for a description of all permanent magnetic material characteristics of this package.
</p>
</html>"));
      record Sm2Co17 =
        Modelica_Magnetic.Material.HardMagnetic.PermanentMagnetData (
          H_cBRef=750000,
          B_rRef=1.02,
          T_refCelsius=20,
          alpha_Br=-0.0003) "Sm2Co17 sintered, exemplary values"
                                      annotation (Documentation(info="<html>
<p>
Please refer to the description of  the enclosing package <a href=\"HardMagnetic\">HardMagnetic</a> for a description of all permanent magnetic material characteristics of this package.
</p>
</html>"));
      record SmCo5 =
        Modelica_Magnetic.Material.HardMagnetic.PermanentMagnetData (
          H_cBRef=720000,
          B_rRef=0.95,
          T_refCelsius=20,
          alpha_Br=-0.0004) "SmCo5 sintered, exemplary values"
                                      annotation (Documentation(info="<html>
<p>
Please refer to the description of  the enclosing package <a href=\"HardMagnetic\">HardMagnetic</a> for a description of all permanent magnetic material characteristics of this package.
</p>
</html>"));
      record PlasticNdFeB =
        Modelica_Magnetic.Material.HardMagnetic.PermanentMagnetData (
          H_cBRef=400000,
          B_rRef=0.6,
          T_refCelsius=20,
          alpha_Br=-0.001) "Plastic-bonded NdFeB, exemplary values"
                                      annotation (Documentation(info="<html>
<p>
Please refer to the description of  the enclosing package <a href=\"HardMagnetic\">HardMagnetic</a> for a description of all permanent magnetic material characteristics of this package.
</p>
</html>"));
      record PlasticSmCo =
        Modelica_Magnetic.Material.HardMagnetic.PermanentMagnetData (
          H_cBRef=385000,
          B_rRef=0.57,
          T_refCelsius=20,
          alpha_Br=-0.0004) "Plastic-bonded Sm-Co, exemplary values"
                                      annotation (Documentation(info="<html>
<p>
Please refer to the description of  the enclosing package <a href=\"HardMagnetic\">HardMagnetic</a> for a description of all permanent magnetic material characteristics of this package.
</p>
</html>"));
      record HardFerrite =
        Modelica_Magnetic.Material.HardMagnetic.PermanentMagnetData (
          H_cBRef=170000,
          B_rRef=0.38,
          T_refCelsius=20,
          alpha_Br=-0.002) "Hard ferrite sintered, exemplary values"
                                      annotation (Documentation(info="<html>
<p>
Please refer to the description of  the enclosing package <a href=\"HardMagnetic\">HardMagnetic</a> for a description of all permanent magnetic material characteristics of this package.
</p>
</html>"));
      record PlasticHardFerrite =
        Modelica_Magnetic.Material.HardMagnetic.PermanentMagnetData (
          H_cBRef=130000,
          B_rRef=0.21,
          T_refCelsius=20,
          alpha_Br=-0.002) "Plastic-bonded hard ferrite, exemplary values"
                                      annotation (Documentation(info="<html>
<p>
Please refer to the description of  the enclosing package <a href=\"HardMagnetic\">HardMagnetic</a> for a description of all permanent magnetic material characteristics of this package.
</p>
</html>"));
  public
      record PermanentMagnetData "Record for permanent magnetic material data"
        extends Modelica.Icons.Record;

        parameter SI.MagneticFieldStrength H_cBRef = 1
        "Coercivity at reference temperature";
        parameter SI.MagneticFluxDensity B_rRef = 1
        "Remanence at reference temperature";
        parameter SI.CelsiusTemperature T_refCelsius = 20
        "Reference temperature";
        parameter Utilities.TemperatureCoefficient alpha_Br = 0
        "Temperature coefficient of remanence";

        annotation (Documentation(info="<html>
<p>
Please refer to the description of  the enclosing package <a href=\"HardMagnetic\">HardMagnetic</a> for a description of all permanent magnetic material characteristics of this package.
</p>
</html>"));
      end PermanentMagnetData;
      annotation (Documentation(info="<html>
<p>
Typical values for remanence, coercivity and the temperature coefficient of remanence are provided for the common permanent magnetic materials illustrated below.
<dl>
<dd>
<p><IMG SRC=\"../Images/Magnetic/Material/HardMagnetic/HardMagneticMaterials.png\" ALT=\"Demagnetization characteristics of included permanent magnetic materials\"></p>
</dd>
</dl>
<p>
Linear demagnetization curves are modelled. The characteristic, temperature-dependent \"knee\" of many permanent magnetic materials is not considered, since proper design of permanent magnetic circuits should avoid operation of permanent magnets \"below\" that point due to partial demagnetization. As a result, the temperature coefficient of coercivity is not considered. Only the temperature coefficient of remanence alpha_Br is accounted for, since it describes the dependende of the demagnetization curve on the temperature sufficiently for the region \"above the knee-point\".
</p>
<p>
Additional user-specific materials can be defined as needed.
</p>
</html>"));
    end HardMagnetic;
  end Material;


  package Sensors "Sensors to measure variables in magnetic networks"

  model MagneticPotentialDifferenceSensor
    "Sensor to measure magnetic potential difference"
    extends Modelica.Icons.RotationalSensor;

    Modelica.Blocks.Interfaces.RealOutput V_mag(
        redeclare type SignalType = SI.MagneticPotentialDifference)
      "Magnetic potential difference between ports p and n as output signal"
       annotation (extent=[-10, -90; 10, -110],
        rotation=90);

    annotation (
      Coordsys(
        extent=[-100, -100; 100, 100],
        grid=[1, 1],
        component=[20, 20]),
      Window(
        x=0.28,
        y=0.29,
        width=0.6,
        height=0.6),
      Icon(Text(
          extent=[-52,1; 48,-57],
          string="V_mag",
          style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=0,
            rgbfillColor={0,0,0})),
        Line(points=[-70,0; -90,0], style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=0,
            rgbfillColor={0,0,0})),
        Line(points=[70,0; 90,0], style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=0,
            rgbfillColor={0,0,0})),
        Line(points=[0,-90; 0,-70]),
        Text(extent=[-140,120; 140,80],   string="%name")),
      Diagram(
        Line(points=[-70,0; -90,0],   style(color=0)),
        Line(points=[70,0; 90,0],   style(color=0)),
        Line(points=[0,-90; 0,-70])));
    Interfaces.PositiveMagneticPort p annotation (extent=[-110,-10; -90,10]);
    Interfaces.NegativeMagneticPort n annotation (extent=[90,-10; 110,10]);
  equation
    p.Phi = 0;
    n.Phi = 0;
    V_mag = p.V_mag - n.V_mag;
  end MagneticPotentialDifferenceSensor;

  model MagneticFluxSensor "Sensor to measure magnetic flux"
    extends Modelica.Icons.RotationalSensor;

    Modelica.Blocks.Interfaces.RealOutput Phi(
        redeclare type SignalType = SI.MagneticFlux)
      "Magnetic flux from port p to port n as output signal"
       annotation (extent=[-10, -90; 10, -110],
        rotation=90);
    annotation (
      Coordsys(
        extent=[-100, -100; 100, 100],
        grid=[1, 1],
        component=[20, 20]),
      Window(
        x=0.23,
        y=0.07,
        width=0.6,
        height=0.6),
      Icon(Text(
          extent=[-29, -11; 30, -70],
          style(color=0),
          string="Phi"),
        Line(points=[-70,0; -90,0],   style(color=0)),
        Text(extent=[-140,120; 140,80],   string="%name"),
        Line(points=[70,0; 90,0],   style(color=0)),
        Line(points=[0,-90; 0,-70])),
      Diagram(
        Line(points=[-70,0; -90,0],   style(color=0)),
        Line(points=[70,0; 90,0],   style(color=0)),
        Line(points=[0,-90; 0,-70])));
    Interfaces.PositiveMagneticPort p annotation (extent=[-110,-10; -90,10]);
    Interfaces.NegativeMagneticPort n annotation (extent=[90,-10; 110,10]);
  equation
    p.V_mag = n.V_mag;
    Phi = p.Phi;
    0 = p.Phi + n.Phi;
  end MagneticFluxSensor;
    annotation (Documentation(info="<html>
<p>
For analysis of magnetic networks, only magnetic potential differences and magnetic flux are variables of interest. For that reason, a magnetic potential sensor is not provided.
</p>
</html>"));
  end Sensors;


  package Utilities "Miscellaneous model components for modelling of actuators"
    extends Modelica.Icons.Library;

    type TemperatureCoefficient = Real (final quantity="TemperatureCoefficient", final unit
        =                                                                                   "1/K");

    class TranslatoryStopper
    "Translatory, nonlinear spring/damper combination with free travel for modelling of impact, e.g. at a stopper"

      extends Modelica.Mechanics.Translational.Interfaces.Compliant;
      parameter Real c(
        final unit="N/m",
        final min=0) = 1e6 "Spring stiffness between impact partners";
      parameter Real d(
        final unit="N/ (m/s)",
        final min=0) = 2e2 "Damping coefficient between impact partners";
      SI.Velocity v_rel "Relative velocity between impact partners";
      Boolean Contact "False, if s_rel > 0";
      parameter SI.Length s_n=1e-5
      "Normalized penetration depth of impact partners";
      parameter Real n=1 "Exponent for force function (n >= 1)";
      SI.Velocity delta_v "Relative velocity";
      annotation (
        Coordsys(
          extent=[-100, -100; 100, 100],
          grid=[2, 2],
          component=[20, 20]),
        Diagram,
        Icon(
          Rectangle(extent=[10,60; 60,-60], style(
              pattern=0,
              fillColor=7,
              rgbfillColor={255,255,255})),
          Text(extent=[-140, 120; 140, 60], string="%name"),
          Polygon(points=[50,-90; 20,-80; 20,-100; 50,-90],     style(color=10,
                fillColor=10)),
          Line(points=[-60,-90; 20,-90],   style(color=0, fillColor=10)),
          Line(points=[10,-60; 10,60], style(color=0, rgbcolor={0,0,0})),
          Line(points=[-60,30; -38,30], style(color=0)),
          Line(points=[-60,-30; -50,-30],  style(color=0)),
          Line(points=[-50,-30; -46,-20; -38,-40; -30,-20; -22,-40; -14,-20;
                -10,-30],
                      style(color=0)),
          Rectangle(extent=[-38,50; -18,10],
                                           style(
              color=0,
              fillColor=8,
              fillPattern=1)),
          Line(points=[-60,30; -60,-30], style(color=0, rgbcolor={0,0,0})),
          Line(points=[-18,30; 10,30], style(color=0, rgbcolor={0,0,0})),
          Line(points=[10,-30; -10,-30], style(color=0, rgbcolor={0,0,0})),
          Line(points=[60,-60; 60,60], style(color=0, rgbcolor={0,0,0})),
          Line(points=[-60,0; -100,0], style(color=0, rgbcolor={0,0,0})),
          Line(points=[60,0; 100,0], style(color=0, rgbcolor={0,0,0})),
          Line(points=[-50,14; -10,54], style(color=0, rgbcolor={0,0,0})),
          Polygon(points=[-8,56; -16,52; -12,48; -8,56],        style(color=10,
                fillColor=10)),
          Line(points=[-50,-50; -10,-10], style(color=0, rgbcolor={0,0,0})),
          Polygon(points=[-8,-8; -16,-12; -12,-16; -8,-8],      style(color=10,
                fillColor=10)),
          Line(points=[-50,50; -38,50], style(color=0, rgbcolor={0,0,0})),
          Line(points=[-50,10; -38,10], style(color=0, rgbcolor={0,0,0}))),
        Window(
          x=0.04,
          y=0.01,
          width=0.6,
          height=0.83),
        Documentation(info="<html>
<p>
This component is intended for modelling of impact between a stopper and a freely travelling translatory mass. It is similar to <a href=\"Modelica.Mechanics.Translational.ElastoGap\">Modelica.Mechanics.Translational.ElastoGap</a>, but the force during contact of the impact partners is calculated differently. Kinetic energy is dissipated in a damper only during collision of the two impact partners, but not if they move away from each other. Also, as penetration depth increases, the force increases non-linearly.
</p>
</html>"));
    equation

      v_rel = der(s_rel);
      Contact = s_rel < 0;
      delta_v = if (v_rel < 0) then v_rel else 0;
      f = if Contact then (c*s_rel + d*delta_v) * (abs(s_rel)/s_n)^n else 0;
    end TranslatoryStopper;

    model TranslatoryArmature "Mass with free travel between two stoppers"

      parameter SI.Mass m=1 "Armature mass";
      parameter SI.Position x_max=10e-3
      "Position of stopper at maximum armature position";
      parameter SI.Position x_min=0
      "Position of stopper at minimum armature position";
      Modelica.Mechanics.Translational.SlidingMass mass(m=m)
        annotation (extent=[-10,-10; 10,10]);
      Modelica.Mechanics.Translational.Interfaces.Flange_a flange_a
        annotation (extent=[-110,-10; -90,10]);
      Modelica.Mechanics.Translational.Interfaces.Flange_b flange_b
        annotation (extent=[90, -10; 110, 10]);
      Modelica.Mechanics.Translational.Fixed limit_xMin(s0=x_min)
        annotation (extent=[-80,-50; -60,-30]);
      Modelica.Mechanics.Translational.Fixed limit_xMax(s0=x_max)
        annotation (extent=[60,-50; 80,-30]);
      Modelica_Magnetic.Utilities.TranslatoryStopper stopper_xMax
        annotation (extent=[50,-30; 70,-10]);
      Modelica_Magnetic.Utilities.TranslatoryStopper stopper_xMin
        annotation (extent=[-70,-30; -50,-10]);

      annotation (
        Coordsys(
          extent=[-100, -100; 100, 100],
          grid=[2, 2],
          component=[20, 20]),
        Window(
          x=0.21,
          y=0.35,
          width=0.53,
          height=0.73),
        Diagram(
          Line(points=[-80,-60; -80,-80], style(color=9, thickness=2)),
          Line(points=[-80,-76; -40,-76], style(
              color=9,
              thickness=2,
              arrow=1)),
          Text(
            extent=[-84,-82; -76,-92],
            string="0",
            style(
              color=9,
              pattern=2,
              thickness=2)),
          Text(
            extent=[-46,-82; -38,-92],
            style(
              color=9,
              pattern=2,
              thickness=2),
            string="x"),
          Ellipse(extent=[-82,-78; -78,-74], style(
              color=9,
              rgbcolor={175,175,175},
              fillColor=9,
              rgbfillColor={175,175,175}))),
        Icon(
          Rectangle(extent=[-90, 10; 90, -10], style(gradient=2, fillColor=45)),
          Text(
            extent=[-120,140; 120,100],
            string="%name",
            style(color=3, rgbcolor={0,0,255})),
          Rectangle(extent=[-50,60; 50,-60],   style(
              color=45,
              gradient=2,
              fillColor=45)),
        Rectangle(extent=[-80,-20; -88,-80], style(
            color=45,
            rgbcolor={255,128,0},
            fillColor=45,
            rgbfillColor={255,128,0})),
        Rectangle(extent=[88,-20; 80,-80], style(
            color=45,
            rgbcolor={255,128,0},
            fillColor=45,
            rgbfillColor={255,128,0})),
        Rectangle(extent=[-80,80; -88,20], style(
            color=45,
            rgbcolor={255,128,0},
            fillColor=45,
            rgbfillColor={255,128,0})),
        Rectangle(extent=[88,80; 80,20], style(
            color=45,
            rgbcolor={255,128,0},
            fillColor=45,
            rgbfillColor={255,128,0})),
          Text(
            extent=[-100,-100; 100,-140],
            style(color=0, rgbcolor={0,0,0}),
            string="m=%m"),
          Line(points=[-50,-80; 30,-80],   style(color=0, fillColor=10)),
          Polygon(points=[60,-80; 30,-70; 30,-90; 60,-80],      style(color=10,
                fillColor=10))),
        Documentation(info="<html>
<p>
In translatory actuators with limited stroke, the armature with its inertia can travel between two stoppers.
</p>
</html>"));
    equation
      connect(mass.flange_a, stopper_xMin.flange_b)
        annotation (points=[-10,0; -40,0; -40,-20; -50,-20],
                                               style(color=58));
      connect(limit_xMax.flange_b,stopper_xMax. flange_b)
        annotation (points=[70,-40; 70,-20], style(color=58));
      connect(stopper_xMax.flange_a, mass.flange_b)
        annotation (points=[50,-20; 40,-20; 40,0; 10,0],
                                             style(color=58));
      connect(mass.flange_a, flange_a)         annotation (points=[-10,0; -100,0],
                        style(color=58, rgbcolor={0,127,0}));
      connect(limit_xMin.flange_b, stopper_xMin.flange_a) annotation (points=[-70,
            -40; -70,-20], style(
          color=58,
          rgbcolor={0,127,0},
          fillColor=9,
          rgbfillColor={175,175,175},
          fillPattern=1));
      connect(flange_b, mass.flange_b) annotation (points=[100,0; 10,0], style(
          color=58,
          rgbcolor={0,127,0},
          fillColor=9,
          rgbfillColor={175,175,175},
          fillPattern=1));
    end TranslatoryArmature;

        block VariableGain "Gain block with variable gain"

  public
          Modelica.Blocks.Interfaces.RealInput u "Input signal connector"
            annotation (extent=[-140, -20; -100, 20]);
          Modelica.Blocks.Interfaces.RealOutput y "Output signal connector"
            annotation (extent=[100, -10; 120, 10]);
          annotation (
            Coordsys(
              extent=[-100, -100; 100, 100],
              grid=[2, 2],
              component=[20, 20]),
            Window(
              x=0.19,
              y=0.02,
              width=0.59,
              height=0.6),
            Documentation(info="<HTML>
<p>
This component is similar to <a href=\"Modelica.Blocks.Math.Gain\">Modelica.Blocks.Math.Gain</a>. However, gain <i>k</i> is not a parameter here and thus has a higher variability. The output <i>y</i> is computed as product of gain <i>k</i> and input <i>u</i>:
<pre>
    y = k * u;
</pre>
</p>
</HTML>
"),         Icon(
              Polygon(points=[-70,-60; -70,60; 70,0; -70,-60],            style(
                  color=74,
                  rgbcolor={0,0,127},
                  fillColor=7,
                  rgbfillColor={255,255,255})),
              Text(extent=[-150, 140; 150, 100], string="%name"),
          Line(points=[-100,0; -70,0], style(color=3, rgbcolor={0,0,255})),
          Line(points=[70,0; 100,0], style(color=3, rgbcolor={0,0,255})),
          Polygon(points=[-6,16; -14,12; -10,8; -6,16], style(
              color=0,
              rgbcolor={0,0,0},
              fillColor=0,
              rgbfillColor={0,0,0})),
          Line(points=[-44,-22; -6,16], style(
              color=0,
              rgbcolor={0,0,0},
              fillColor=0,
              rgbfillColor={0,0,0}))),
            Diagram(Polygon(points=[-100, -100; -100, 100; 100, 0; -100, -100],
                      style(
                        color=74,
                        rgbcolor={0,0,127},
                        fillColor=7,
                        rgbfillColor={255,255,255}))));
          Modelica.Blocks.Interfaces.RealInput k "Variable gain"
            annotation (extent=[-20,-70; 20,-30], rotation=90);
        equation
          y = k*u;
        end VariableGain;

  end Utilities;


  model MagneticGround "Zero magnetic potential"

    Modelica_Magnetic.Interfaces.MagneticPort p
      annotation (extent=[-10, 110; 10, 90], rotation=-90);
    annotation (
      Coordsys(
        extent=[-100, -100; 100, 100],
        grid=[2, 2],
        component=[20, 20]),
      Documentation(info="<HTML>
<P>
The magnetic potential at the magnetic ground node is zero. Every magnetic network model must contain at least one magnetic ground object.
</P>
</HTML>
"),   Icon(
        Line(points=[-60, 50; 60, 50], style(color=45)),
        Line(points=[-40, 30; 40, 30], style(color=45)),
        Line(points=[-20, 10; 20, 10], style(color=45)),
        Line(points=[0, 90; 0, 50], style(color=45)),
        Text(
          extent=[-100,-40; 100,-2],
          string="%name",
          style(color=3, rgbcolor={0,0,255}))),
      Diagram(
        Line(points=[-60, 50; 60, 50], style(color=45, thickness=2)),
        Line(points=[-40, 30; 40, 30], style(color=45, thickness=2)),
        Line(points=[-20, 10; 20, 10], style(color=45, thickness=2)),
        Line(points=[0,100; 0,50],  style(color=45, thickness=2)),
        Text(
          extent=[-40,-40; 40,20],
          string="p.V_mag = 0",
          style(color=3, rgbcolor={0,0,255}))),
      Window(
        x=0.23,
        y=0.23,
        width=0.59,
        height=0.63));
  equation
    p.V_mag = 0;
  end MagneticGround;
end Modelica_Magnetic;
