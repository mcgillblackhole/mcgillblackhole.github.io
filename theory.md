---
title: Theory
---

<h2>Physics behind the project</h2>

<html>
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width">
  <title>MathJax example</title>
  <script src="https://polyfill.io/v3/polyfill.min.js?features=es6"></script>
  <script id="MathJax-script" async
          src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js">
  </script>
</head>
<body>
<p>
  <p style="text-align:justify">The three-body problem is famous for being a chaotic system without analytical solution. A simple example of three-body systems is
  the description of the Earth-Sun-Jupiter orbits. Usually in this type of situation, physicists ignore the mass of the Earth, since both
  Jupiter and the Sun are more massive. Another three-body system that is usually approximated to a two body system is the Alpha Centauri star system.
  The system is composed of Alpha Centauri A, B and C. However, Alpha Centauri C (also known as Proxima Centauri) mass is only 0.122 solar masses
  (denoted by \(M_\odot\) ), while the masses of Alpha Centauri A and B are 1.1 \(M_\odot\) and 0.907 \(M_\odot\), respectively. Furthermore, Proxima Centauri
  appears to be far away from the first two, not influencing their dynamics significantly. The system is then approximated to a binary sistem of Alpha
  Centauri A and B (See ref [<a href="https://arxiv.org/abs/astro-ph/9609106"> 1</a>] for more details).<br> </p>

 <p style="text-align:justify"> In this project, we do not wish to work with this approximation. That means we need to solve numerically the system of coupled equations in order to extract
  the orbits of the bodies. The equation of motion of each \(i\)-th body in the system is given by the Einstein-Infeld-Hoffmann equation [<a href="https://doi.org/10.2307%2F1968714">2</a>]
  
  \[ \vec{a}_i = \sum_{i \neq j}{G m_j \vec{n}_{ij} \over r_{ij}^2} + {1 \over c^2} \sum_{i \neq j} {G m_j \vec{n}_{ij} \over r_{ij}^2} {\left[v_{i}^2 + 2v_{j}^2
  - 4 \left( \vec{v}_i \cdot \vec{v}_j \right) - {3 \over 2} \left( \vec{n}_{ij} \cdot \vec{v}_j \right)^2 - 4 \sum_{k \neq i} {G m_k \over r_{ik}} - \right.}\]
  \[{\left. - \sum_{j \neq k} {G m_k \over r_{jk}} + {1 \over 2} \left( \left( \vec{r}_j - \vec{r}_i \right) \cdot \vec{a}_j \right) \right]} + {1 \over c^2}
  \sum_{i \neq j} {G m_j \over r_{ij}^2} \left[ \vec{n}_{ij} \cdot \left(4 \vec{v}_i - 3 \vec{v}_j \right) \right] \left( \vec{v}_i - \vec{v}_j \right) + \]
  \[+ {7 \over 2 c^2} \sum_{i \neq j} {G m_j \vec{a}_{j} \over r_{ij}^2} + \mathcal{O}\left( c^{-4} \right)\] </p>
  
 <p style="text-align:justify"> This is a pretty complicated equation, the acceleration of the \(i\)-th body depends on the acceleration of the others \(j\)-th bodies and it
  also depends on the velocities and the positions of all other \(j\)-th bodies, where \(i \neq j\). In our case, the total number of orbiting bodies is 3. The
  \(\vec{n}_{ij}\) vector is the unit vector pointing from the \(i\)-th body to the \(j\)-th body. We then have 3 equations like this one, one for each body.<br>
   <br> </p>
  
<p style="text-align:justify">  In order to solve all three equations, we set initial conditions and the constant parameters (\(G, m_1, m_2, m_3\)). Since this is a relativistic problem, we may
  expect some gravitation waves being emitted from the orbital motion of the three bodies (just like the two black holes merger that provided the first gravitational
  wave detection [<a target="_blank" href="https://doi.org/10.1103/PhysRevLett.116.061102">3</a>]). The emission of gravitational waves close to the emitting source
  is a very complex problem. In our project, we are interested in how we see this gravitational waves from a frame of reference far way from the source. We may think
  of this consideration as equivallent to measuring gravitational waves here on Earth from the other side of the galaxy. Therefore, we may consider only the first
  order of approximation from the multipole expansion using perturbation theory.
  
  When considering the emission of gravitational waves, we are talking about a perturbation on the metric tensor \(g_{\mu \nu}\) that describes space-time
  
  \[G_{\mu \nu} = {8 \pi \over c^4} T_{\mu \nu} \]
  
  where the metric is contained inside \(G_{\mu \nu}\). We first...
  
  \[g_{\mu \nu} = \eta_{\mu \nu} + h_{\mu \nu} \]
  
  Then mathemagic happens and boom </p>
<p>
<body>


<h1>References</h1>

<ol>
  <li> Wiegert, P. A., & Holman, M. (1996). The stability of planets in the Alpha Centauri system. arXiv preprint <a href="https://arxiv.org/abs/astro-ph/9609106"> astro-ph/9609106</a>;</li>
  <li>Einstein, A., Infeld, L., & Hoffmann, B. (1938). The gravitational equations and the problem of motion. Annals of mathematics, 65-100. <a href="https://doi.org/10.2307%2F1968714"> doi:10.2307/1968714</a>;</li>
  <li>Abbott, B. P., Abbott, R., Abbott, T. D., Abernathy, M. R., Acernese, F., Ackley, K., ... & Cavalieri, R. (2016). Observation of gravitational waves from a binary black hole merger. 
    Physical review letters, 116(6), 061102.<a target="_blank" href="https://doi.org/10.1103/PhysRevLett.116.061102">doi:10.1103/PhysRevLett.116.061102.</a></li>
  <li>Li, X., Jing, Y., & Liao, S. (2018). Over a thousand new periodic orbits of a planar three-body system with unequal masses. Publications of the Astronomical
    Society of Japan, 70(4), 64.<a target="_blank" href="https://doi.org/10.1093/pasj/psy057">doi:10.1093/pasj/psy057.</a></li>
  <li>Dmitrašinović, V., Šuvakov, M., & Hudomal, A. (2014). Gravitational waves from periodic three-body systems. Physical review letters, 113(10), 101102.<a target="_blank" href="https://doi.org/10.1103/PhysRevLett.113.101102">doi:10.1103/PhysRevLett.113.101102</a></li>
</ol>
