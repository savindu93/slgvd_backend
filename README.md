# Sri Lanka Genetic Variation Database (SLGVD)

SLGVD is the only primary database for genetic variant 🧬 data on Sri Lankans hosted by the Faculty of Medicine, University of Colombo. 

The 2nd iteration of SLGVD serves as an archive for single nucleotide variants (SNVs), short indels, and germline copy number variants (gCNVs). It also allows users to search for variants using complex queries. Moreover, it provides a platform for users to download ⬇️, upload ⬆️, edit 🖋️ and remove variant data.

The web application was developed as,
- A React application with
- A back-end implemented using Django and
- A SQL database that serves as the primary storage unit for the variation data at SLGVD.

The Django Rest Framework provides a public REST API for users to programmetically access the data in the database. 

The following repo presents the code base behind the backend of the SLGVD.
