# Waning notes

We may want to add in when waning stops. For incubation times we need to draw $q \sim U(n, 0, 1)$. Under a given waning function, $w(t)$, find $\tau = \arg\min\limits_t \{w(t) < q\}$. Modify $1 - w(t)$ to be $1 - w(t)\mathbf{I}(w(t) > q) \equiv 1 - w(t)\mathbf{I}(t > \tau)$ .