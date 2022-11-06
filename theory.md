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
  <p style="text-align:justify">Black Holes have always been objects of interest to everyone, from scientists to lovers of science fiction. Even though these objects were proposed some decadesa go, it was only recently that we were finnaly able to observe them in the Event Horizon Collaboration, and with these new images, the interest in them seems to have been rekindled, our interest included. So, we decided to work trying to simulate the photograph of a Black Hole, something that was already proposed many years ago (as a reference, check [<a href="https://arxiv.org/abs/astro-ph/9609106">1</a>]) and that can reproduce results very similar to the ones captured by the telescopes here on Earth.</p>
<br>
 <p style="text-align:justify"> In this project, we consider only Schwarszchild Black Holes, with metric given by
   
   $$ds^2 = -\left(1-\frac{2M}{r}\right)dt^{2} + \left(1-\frac{2M}{r}\right)^{-1} + r^{2}(d\theta^{2} + sin^{2}(\theta)d\phi^{2}).$$
   
   The first thing that had to be done was to find a way and understand how to solve the equations of motion for a massless particle (photon), since this would be the base for all the rest of our work. So, from the metric we can get the following equations (we could also use the Lagrangean to find the equations of motion), using \(\theta = \pi/2\) and having \(b\) as an impact parameter

   $$\left\{\frac{1}{r^{2}}\left(\frac{dr}{d\phi}\right)^{2}\right\} + \frac{1}{r^{2}}\left(1-\frac{2M}{r}\right) = \frac{1}{b^{2}},$$
   
and with \(u = 1/r\)</p>
   
   $$\left(\frac{du}{d\phi}\right)^{2} = 2Mu^{3} - u^{2} + \frac{1}{b^{2}}.$$
   
  <p style="text-align:justify">Solving this equation we can find the lightlike geodesics we are interested in, as shown in the figure bellow. Following our reference, we can also find some expressions useful for finding level curves (contours) of the accretion disk surrounding our Black Hole. To do this, we need to perform some quite lenghty calculations and to define new variables \(b = P^{3} / (P-2M)\), \(Q^{2} = (P-2M)(P+6M)\), that relates our impact parameter \(b\) to the periastron distance \(P\). With this, and solving some elliptical integrals, we can get
  
   $$\frac{1}{r} = -\frac{Q - P + 2M}{4MP} + \frac{Q - P + 6M}{4MP} sn^{2}\left\{\frac{\gamma}{2}\sqrt{\frac{Q}{P}} + F(\xi_{\infty}, k)\right\},$$
    
   and use this expression to find contours as the following.</p>

  

   $$F_{s} = \frac{3M\dot{M}}{8\pi}\frac{1}{(r^* - 3)r^{*5/2}}\left[ \sqrt{r^{*}} - \sqrt{6} + \frac{\sqrt{3}}{3}log\left\{ \frac{(\sqrt{r^{*}} + \sqrt{3})(\sqrt{6} - \sqrt{3})}{(\sqrt{r^{*}} - \sqrt{3})(\sqrt{6} + \sqrt{3})} \right\} \right]$$

   $$1+z = (1-3\frac{M}{r})^{-1/2}(1+\left(\frac{M}{r^{3}}\right)^{1/2}b\sin(\theta)\sin(\alpha)$$
  
  
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
