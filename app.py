import streamlit as st

from textbook import CHAPTERS

st.set_page_config(page_title="线性代数", layout="wide")

st.sidebar.title("线性代数")
titles = [c["title"] for c in CHAPTERS]
choice = st.sidebar.radio("章节", titles, label_visibility="collapsed")
chapter = CHAPTERS[titles.index(choice)]

st.title(chapter["title"])
st.markdown(chapter["theory_md"])

st.header("习题")
for i, exercise in enumerate(chapter["exercises"], 1):
    st.markdown(f"**习题 {i}**")
    st.markdown(exercise["problem_md"])
    with st.expander("解答"):
        st.markdown(exercise["solution_md"])
