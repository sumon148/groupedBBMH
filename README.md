
# groupedBBMH

The **groupedBBMH** R package implements a **nested beta-binomial
hierarchical model** designed for **group-testing data in biosecurity
surveillance**, where the primary goal is to efficiently detect rare
contamination events. In routine inspections of imported agricultural
consignments, only group-level binary outcomes—positive or negative—are
observed, without information on the exact number of contaminated items.
This package extends traditional beta-binomial modeling to incorporate
**imperfect test sensitivity and specificity**, offering both **exact
inference** using a Metropolis-Hastings (MH) algorithm and **fast
approximate estimation** via a beta-binomial approximation. It supports
estimation of contamination levels, quantifies the risk of undetected
contamination, and aids in risk-based decision-making through
model-based simulation.

To illustrate the method, consider a biosecurity scenario involving the
importation of frozen seafood into a country. Suppose there are $D$
consignments (or batches), each consisting of approximately $B$ groups,
with each group containing $m$ items — so that the total number of items
per batch is $N = mB$.

From each batch, a simple random sample of $b$ groups is selected for
group testing. These tests, which screen for pathogens or contaminants,
are typically imperfect. We assume a constant sensitivity $\Delta$ — the
probability that a test correctly detects infection when it is present —
and a constant specificity $\Lambda$ — the probability that a test
correctly identifies no infection when it is absent.

For each batch, the only observed outcome is the number of positive
groups among the $b$ tested.

Let $X_{ij\ell}$ indicate contamination status (1 = contaminated, 0 =
not) for unit $\ell$ in group $j$ of consignment $i$, where
$j = 1, \ldots, b$ and $\ell = 1, \ldots, m$. The number of contaminated
units in group $j$ of consignment $i$ is
$X_{ij} = \sum_{\ell=1}^{m} X_{ij\ell}$, taking values
$0, 1, \ldots, m$.

The values $X_{ij}$ are not directly observed. Instead, we observe
$Y_{ij} = I(X_{ij} > 0)$, indicating whether any contamination is
present in group $j$. We then define:

- $T_{Xi} = \sum_{j=1}^B X_{ij}$: total contaminated units in
  consignment $i$,
- $t_{xi} = \sum_{j=1}^{b} X_{ij}$: contaminated units in the sampled
  groups,
- $t_{yi} = \sum_{j=1}^{b} Y_{ij}$: number of sampled groups with
  contamination, ranging from 0 to $b$.

Let $\tilde{Y}_{ij}$ denote the observed test result for group $j$ in
consignment $i$, with 1 indicating a positive result (i.e., apparent
contamination) and 0 denoting a negative result (apparently
contamination-free). The probability of a positive test results is then

$$
\Pr \left( \tilde{Y}_{ij} = 1 \vert Y_{ij} \right) = \Delta Y_{ij} + \left( 1-\Lambda \right) \left(1-Y_{ij}\right). 
$$

Thus, the observed data for each consignment consists of
$\tilde{t}_{yi}=\sum_{j=1}^{b} \tilde{Y}_{ij}$.

A key quantity of interest is , defined as:

$$
L_i = T_{Xi} \cdot \mathbb{I}(\tilde{t}_{yi} = 0),
$$

representing the number of contaminated items entering undetected when
no sampled group tests positive. We focus on:

- the **expected leakage** $\mathbb{E}(L_i)$,
- the **probability of leakage** $\Pr(L_i > 0)$,
- and posterior inference on $T_{Xi}$ when $\tilde{t}_{yi} = 0$.

Now let $p_i$ denote the true prevalence of contamination in batch $i$,
and define $\phi_{i} = 1 - (1 - p_i)^m$ as the probability that a group
of $m$ items is contaminated. Assuming items are randomly distributed
among bags, we have:

$$
X_{ij} \mid p_i \overset{\text{i.i.d.}}{\sim} \text{Bin}(m, p_i), \quad p_i \overset{\text{i.i.d.}}{\sim} \text{Beta}(\alpha, \beta).
$$

The group-level test outcome $\tilde{Y}_{ij}$ follows:

$$
\tilde{Y}_{ij} \mid p_i \sim \text{Bernoulli}(\tilde{\phi}_i), 
$$

and the observed number of positive groups in the sample follows:

$$
\tilde{t}_{yi} \mid p_i \sim \text{Binomial}(b, \tilde{\phi}_i).
$$

where the effective probability of a positive test is:

$$
\tilde{\phi} = \Delta \phi_i + (1 - \Lambda)(1 - \phi_i) \quad \phi_i = 1 - (1 - p_i)^m.
$$

Under perfect specificity ($\Lambda = 1$), this simplifies to:

$$
\tilde{\phi}_i = \Delta \phi_i.
$$

When $\beta \gg \alpha$, the contamination prevalence $\beta p_i$ is
approximately Gamma distributed, leading to:

$$
\tilde{\phi}_i \sim (1 - \Lambda) + \text{Beta} \left( \alpha, \frac{\beta}{m(\Delta + \Lambda - 1)} \right).
$$

Two important special cases:

- **Perfect testing** ($\Delta = 1$, $\Lambda = 1$):
  $\tilde{\phi}_i \sim \text{Beta}(\alpha, \beta / m)$,
- **Perfect specificity** ($\Lambda = 1$):
  $\tilde{\phi}_i \sim \text{Beta}(\alpha, \frac{\beta}{m \Delta})$.

In either case, $\tilde{t}_{yi}$ approximately follows a Beta-Binomial
distribution:

$$
\tilde{t}_{yi} \sim \textbf{Beta-Binomial}\ \Bigg(b, \alpha, \frac{\beta}{m \Delta}\Bigg).
$$

#### Threshold-based Risk Estimation

In the context of the frozen seafood biosecurity study, suppose a
regulatory threshold is set such that contamination levels below a
certain prevalence cut-off are considered acceptable for import. Suppose
a regulatory cut-off $k$ is introduced, such that contamination
prevalence below $k$ is considered acceptable. The effective prevalence
becomes:

$$
p_i = p_i \cdot \mathbb{I}(p_i > k), \quad p_i \sim \text{Beta}(\alpha, \beta).
$$

Using this truncated model, the **probability of leakage** becomes:

$$
\Pr[L_i > 0] = \frac{B_{(k_1,1)}(\alpha, \frac{\beta}{m\Delta} + b)}{B(\alpha, \frac{\beta}{m\Delta})} - \frac{B_{(k,1)}(\alpha, \beta + MB)}{B(\alpha, \beta)},
$$

with $k_1 = \Delta(1 - (1 - k)^m)$. Under perfect testing:

$$
\Pr[L_i > 0] = \frac{B_{(k,1)}(\alpha, \beta + mb)}{B(\alpha, \beta)} - \frac{B_{(k,1)}(\alpha, \beta + MB)}{B(\alpha, \beta)}.
$$

When $k = 0$, this reduces to the standard beta-binomial case.

The **expected leakage** under threshold $k$ is:

$$
\mathbb{E}(L_i) = (B - b) M \cdot \mathbb{E}\left[ (1 - p_i)^{bm\Delta} \cdot p_i \cdot \mathbb{I}(p_i > k) \right],
$$

which simplifies to:

$$
\mathbb{E}(L_i) = (B - b) M \cdot \frac{B_{(k,1)}(\alpha + 1, \beta + bm\Delta)}{B(\alpha, \beta)}.
$$

#### Model-based Simulation to Assess Leakage Risk

Using the posterior estimates of Beta distribution parameters
$(\alpha, \beta)$, we simulate contamination and testing outcomes for
$D$ consignments of frozen seafood items. Each consignment consists of
$B$ bags, with each bag containing $m$ items. For inspection, a simple
random sample of $b$ bags is selected from each consignment and tested
using PCR-type assays.

For each consignment $i$, the contamination prevalence $p_i$ is
simulated from $\text{Beta}(\alpha, \beta)$ and is used to simulate the
number of contaminated items per bag as

$$
X_{ij} \sim \text{Binomial}(M, p_i)
$$

which provides the true contamination status of bag $j$ as
$Y_{ij} = \mathbb{I}(X_{ij} > 0)$ (i.e., contains at least one
contaminated item). The total number of contamination bags in the
consignment $i$ is

$$
T_{yi} = \sum_{j=1}^{B} Y_{ij}
$$

and the total number of contamination bags in the sample of size $b$ is

$$
{t}_{yi} = \sum_{j=1}^{b} {Y}_{ij}.
$$

To account for imperfect testing, we simulate observed test results
$\tilde{t}_{yi}$ by applying the assumed test sensitivity and
specificity to the sampled contamination statuses ${Y}_{ij}$. The
observed number of positive test results among the $b$ sampled bags is
then

$$
\tilde{t}_{yi} = \sum_{j=1}^{b} \tilde{Y}_{ij}.
$$

Test performance is evaluated by comparing true contamination status of
each consignment ($T_{yi}$) with its observed test outcome
$\tilde{t}_{yi}$. Specifically, consignments are classified as:

- **True Positive (TP):** $T_{yi} > 0$ and $\tilde{t}_{yi} > 0$
- **True Negative (TN):** $T_{yi} = 0$ and $\tilde{t}_{yi} = 0$
- **False Negative (FN):** $T_{yi} > 0$ and $\tilde{t}_{yi} = 0$

False negatives lead to **leakage**, meaning contaminated consignments
pass undetected into the country.

This simulation framework allows estimation of key biosecurity risk
metrics:

- **Probability of leakage**:
  $\Pr(T_{Xi} > 0 \ \text{and} \ \tilde{t}_{yi} = 0)$,
- **Expected leakage**: the average number of contaminated items
  undetected per consignment.

These metrics help to quantify the effectiveness of group testing
strategies and inform decision-making in biosecurity risk management.

#### Case Study: Frozen Seafood Importation in Australia

A case study, drawn from a real biosecurity context, serves as the
primary motivation for our modeling approach. The data are included in
this package for illustration and reproducibility.

Each year, approximately 800 consignments (batches) of frozen seafood
are imported into Australia. Each batch contains roughly $B = 8000$
bags, with each bag holding $M = 40$ items, yielding a total of $N = MB$
items per batch.

A simple random sample of $b = 13$ bags is selected from each batch for
testing. From each selected bag, a group of $m = 5$ items (out of
$M = 40$) is randomly chosen, resulting in $n = bm = 65$ sampled items
per batch. Each group of $m$ items is tested using a PCR assay, which is
assumed to have a sensitivity of approximately $\Delta = 0.80$ and
near-perfect specificity, $\Lambda = 1.0$.

The only observation recorded per batch is the number of groups testing
positive, out of the $b = 13$ groups tested.

These data motivate the development of models for estimating
contamination prevalence, quantifying the risk of undetected
contamination (leakage), and evaluating testing strategies under
imperfect diagnostic accuracy.

## Installation

You can install the development version of `groupedBBMH` like so:

``` r
install.packages("devtools")
devtools::install_github("sumon148/groupedBBMH")
```
