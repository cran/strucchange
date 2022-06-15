\name{scPublications}
\alias{scPublications}
\title{Structural Change Publications}

\usage{data("scPublications")}

\description{
Bibliographic information about papers related to structural change and
changepoints published in 27 different econometrics and statistics 
journals.
}
\format{
A data frame containing information on 835 structural change papers in 9 variables.
 \describe{
    \item{author}{character. Author(s) of the paper.}
    \item{title}{character. Title of the paper.}
    \item{journal}{factor. In which journal was the paper published?}
    \item{year}{numeric. Year of publication.}
    \item{volume}{numeric. Journal volume.}
    \item{issue}{character. Issue within the journal volume.}
    \item{bpage}{numeric. Page on which the paper begins.}
    \item{epage}{numeric. Page on which the paper ends.}
    \item{type}{factor. Is the journal an econometrics or statistics journal?}
}
}

\details{
The data set \code{scPublications} includes
bibliographic information about publications related to structural change and obtained
from the \sQuote{ISI Web of Science}. The query was based on the \sQuote{Science Citation Index Expanded}
and \sQuote{Social Sciences Citation Index} (for the full range of years available: 1900-2006 and
1956-2006, respectively). The \sQuote{Source Title} was restricted to the 27 journals
in the data frame and the \sQuote{Topic} to be one of the following:
structural change, structural break, structural stability, structural instability,
parameter instability, parameter stability, parameter constancy, change point,
changepoint, change-point, breakpoint, break-point, break point, CUSUM, MOSUM.
Additionally, the famous CUSUM paper of Brown, Durbin and Evans (1975) was added
manually to \code{scPublications} (because it did not match the query above).
}

\source{ISI Web of Science at \url{https://www.webofknowledge.com/}.
Queried by James Bullard.
}

\examples{
## construct time series:
## number of sc publications in econometrics/statistics
data("scPublications")

## select years from 1987 and
## `most important' journals
pub <- scPublications
pub <- subset(pub, year > 1986)
tab1 <- table(pub$journal)
nam1 <- names(tab1)[as.vector(tab1) > 9] ## at least 10 papers
tab2 <- sapply(levels(pub$journal), function(x) min(subset(pub, journal == x)$year))
nam2 <- names(tab2)[as.vector(tab2) < 1991] ## started at least in 1990
nam <- nam1[nam1 \%in\% nam2]
pub <- subset(pub, as.character(journal) \%in\% nam)
pub$journal <- factor(pub$journal)
pub_data <- pub

## generate time series
pub <- with(pub, tapply(type, year, table))
pub <- zoo(t(sapply(pub, cbind)), 1987:2006)
colnames(pub) <- levels(pub_data$type)

## visualize
plot(pub, ylim = c(0, 35))
}

\keyword{datasets}
